//! Bit-exact reimplementation of `np.random.default_rng(seed)` (SeedSequence →
//! PCG64) and of the column permutation WitChi performs with
//! `np.apply_along_axis(rng.permutation, 0, alignment_array)`.
//!
//! References: numpy/random/bit_generator.pyx (SeedSequence),
//! numpy/random/src/pcg64 (PCG64 XSL-RR, 128-bit LCG), numpy's
//! `random_interval` (masked rejection sampling, 32-bit draws for bounds
//! below 2^32 using the buffered half of a 64-bit draw) and
//! `Generator.shuffle` (Fisher–Yates from the end).

const INIT_A: u32 = 0x43b0_d7e5;
const MULT_A: u32 = 0x931e_8875;
const INIT_B: u32 = 0x8b51_f9dd;
const MULT_B: u32 = 0x58f3_8ded;
const MIX_MULT_L: u32 = 0xca01_f9dd;
const MIX_MULT_R: u32 = 0x4973_f715;
const XSHIFT: u32 = 16;
const POOL_SIZE: usize = 4;

/// `SeedSequence(entropy).generate_state(n_words, np.uint64)` for a
/// non-negative integer seed.
pub fn seed_sequence_state(entropy: u128, n_words64: usize) -> Vec<u64> {
    // Integer entropy is split into little-endian 32-bit words (at least one).
    let mut ent = Vec::new();
    let mut e = entropy;
    loop {
        ent.push(e as u32);
        e >>= 32;
        if e == 0 {
            break;
        }
    }
    let mut hash_const = INIT_A;
    let mut hashmix = |value: u32| -> u32 {
        let mut v = value ^ hash_const;
        hash_const = hash_const.wrapping_mul(MULT_A);
        v = v.wrapping_mul(hash_const);
        v ^ (v >> XSHIFT)
    };
    let mix = |x: u32, y: u32| -> u32 {
        let r = MIX_MULT_L.wrapping_mul(x).wrapping_sub(MIX_MULT_R.wrapping_mul(y));
        r ^ (r >> XSHIFT)
    };
    let mut pool = [0u32; POOL_SIZE];
    for (i, p) in pool.iter_mut().enumerate() {
        *p = hashmix(ent.get(i).copied().unwrap_or(0));
    }
    for i_src in 0..POOL_SIZE {
        for i_dst in 0..POOL_SIZE {
            if i_src != i_dst {
                let h = hashmix(pool[i_src]);
                pool[i_dst] = mix(pool[i_dst], h);
            }
        }
    }
    for &word in ent.iter().skip(POOL_SIZE) {
        for i_dst in 0..POOL_SIZE {
            let h = hashmix(word);
            pool[i_dst] = mix(pool[i_dst], h);
        }
    }
    let mut out32 = Vec::with_capacity(n_words64 * 2);
    let mut hash_const = INIT_B;
    for i in 0..n_words64 * 2 {
        let mut d = pool[i % POOL_SIZE] ^ hash_const;
        hash_const = hash_const.wrapping_mul(MULT_B);
        d = d.wrapping_mul(hash_const);
        out32.push(d ^ (d >> XSHIFT));
    }
    out32
        .chunks(2)
        .map(|c| c[0] as u64 | ((c[1] as u64) << 32))
        .collect()
}

/// numpy's PCG64 bit generator, including the buffered 32-bit output.
#[derive(Clone, Debug)]
pub struct Pcg64 {
    state: u128,
    inc: u128,
    has_uint32: bool,
    uinteger: u32,
}

const PCG_MULT: u128 = (2549297995355413924u128 << 64) | 4865540595714422341u128;

impl Pcg64 {
    /// `np.random.default_rng(seed)` / `np.random.PCG64(seed)`.
    pub fn from_seed(seed: u64) -> Self {
        let v = seed_sequence_state(seed as u128, 4);
        let initstate = ((v[0] as u128) << 64) | v[1] as u128;
        let initseq = ((v[2] as u128) << 64) | v[3] as u128;
        let mut rng = Pcg64 { state: 0, inc: (initseq << 1) | 1, has_uint32: false, uinteger: 0 };
        rng.step();
        rng.state = rng.state.wrapping_add(initstate);
        rng.step();
        rng
    }

    /// Advance the LCG by `delta` steps in O(log delta) (`pcg_advance_lcg_128`,
    /// numpy's `PCG64.advance`, without touching the buffered half).
    pub fn advance(&mut self, delta: u128) {
        let (mut cur_mult, mut cur_plus) = (PCG_MULT, self.inc);
        let (mut acc_mult, mut acc_plus) = (1u128, 0u128);
        let mut d = delta;
        while d > 0 {
            if d & 1 == 1 {
                acc_mult = acc_mult.wrapping_mul(cur_mult);
                acc_plus = acc_plus.wrapping_mul(cur_mult).wrapping_add(cur_plus);
            }
            cur_plus = cur_mult.wrapping_add(1).wrapping_mul(cur_plus);
            cur_mult = cur_mult.wrapping_mul(cur_mult);
            d >>= 1;
        }
        self.state = acc_mult.wrapping_mul(self.state).wrapping_add(acc_plus);
    }

    #[inline]
    fn step(&mut self) {
        self.state = self.state.wrapping_mul(PCG_MULT).wrapping_add(self.inc);
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.step();
        let s = self.state;
        let x = ((s >> 64) as u64) ^ (s as u64);
        x.rotate_right((s >> 122) as u32)
    }

    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        if self.has_uint32 {
            self.has_uint32 = false;
            return self.uinteger;
        }
        let n = self.next_u64();
        self.has_uint32 = true;
        self.uinteger = (n >> 32) as u32;
        n as u32
    }

    /// numpy `random_interval(bitgen, max)`: uniform integer in [0, max].
    #[inline]
    pub fn random_interval(&mut self, max: u64) -> u64 {
        if max == 0 {
            return 0;
        }
        let mask = u64::MAX >> max.leading_zeros();
        if max <= u32::MAX as u64 {
            loop {
                let v = (self.next_u32() as u64) & mask;
                if v <= max {
                    return v;
                }
            }
        }
        loop {
            let v = self.next_u64() & mask;
            if v <= max {
                return v;
            }
        }
    }

    /// `Generator.shuffle` of a 1-D buffer.
    #[inline]
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.random_interval(i as u64) as usize;
            items.swap(i, j);
        }
    }
}

/// numpy's stream of 32-bit draws from a PCG64 (`next_uint32`: the low then the
/// high half of each 64-bit output), generated in batches, with a branch-free
/// `Generator.shuffle` on top.
///
/// Use it in place of [`Pcg64::shuffle`] when a generator does nothing but
/// shuffle: it consumes exactly the same draws in the same order, so it gives
/// the same permutations and leaves the stream at the same point.
pub struct U32Stream {
    rng: Pcg64,
    buf: Box<[u32; Self::BATCH]>,
    pos: usize,
    len: usize,
    /// Draws put into `buf` so far, counted from the stream start.
    generated: u64,
    /// Scratch for the swap targets of [`shuffle`](Self::shuffle).
    targets: Vec<u32>,
}

impl U32Stream {
    const BATCH: usize = 512;

    pub fn new(mut rng: Pcg64) -> Self {
        let mut stream = U32Stream {
            rng: rng.clone(),
            buf: Box::new([0; Self::BATCH]),
            pos: 0,
            len: 0,
            generated: 0,
            targets: Vec::new(),
        };
        // A pending buffered half comes first, as in next_uint32.
        if rng.has_uint32 {
            rng.has_uint32 = false;
            stream.buf[0] = rng.uinteger;
            stream.len = 1;
            stream.generated = 1;
            stream.rng = rng;
        }
        stream
    }

    /// `np.random.default_rng(seed)`, used only for 32-bit draws.
    pub fn from_seed(seed: u64) -> Self {
        U32Stream::new(Pcg64::from_seed(seed))
    }

    /// The stream of `default_rng(seed)` after its first `position` 32-bit
    /// draws, reached by jumping the LCG ahead instead of drawing. Draw
    /// 2k is the low half of 64-bit output k and draw 2k + 1 its high half.
    pub fn at(seed: u64, position: u64) -> Self {
        let mut rng = Pcg64::from_seed(seed);
        rng.advance((position / 2) as u128);
        let mut stream = U32Stream::new(rng);
        stream.generated = position - position % 2;
        if position % 2 == 1 {
            // Mid-output: the high half of output k is the next draw.
            let x = stream.rng.next_u64();
            stream.buf[0] = (x >> 32) as u32;
            stream.len = 1;
            stream.generated += 2;
        }
        stream
    }

    /// An independent copy of the stream at its current position.
    pub fn clone_stream(&self) -> Self {
        U32Stream {
            rng: self.rng.clone(),
            buf: self.buf.clone(),
            pos: self.pos,
            len: self.len,
            generated: self.generated,
            targets: Vec::new(),
        }
    }

    /// Number of 32-bit draws consumed since the stream start (meaningful
    /// for streams made by [`from_seed`](Self::from_seed) or [`at`](Self::at)).
    pub fn position(&self) -> u64 {
        self.generated - (self.len - self.pos) as u64
    }

    #[inline]
    fn refill(&mut self) {
        for k in 0..Self::BATCH / 2 {
            let x = self.rng.next_u64();
            self.buf[2 * k] = x as u32;
            self.buf[2 * k + 1] = (x >> 32) as u32;
        }
        self.pos = 0;
        self.len = Self::BATCH;
        self.generated += Self::BATCH as u64;
    }

    /// The next 32-bit draw (`next_uint32`).
    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        if self.pos == self.len {
            self.refill();
        }
        let v = self.buf[self.pos];
        self.pos += 1;
        v
    }

    /// `Generator.shuffle` of a 1-D buffer: Fisher–Yates from the end with
    /// `random_interval(i)` (mask rejection on 32-bit draws).
    ///
    /// Two phases. First the swap targets: for each position i = n-1 .. 1 the
    /// first draw whose masked value is <= i, which is exactly what
    /// `random_interval(i)` returns. Within a power-of-two range the mask is
    /// constant, so each draw costs a mask, a compare and a subtract, with no
    /// unpredictable branch (a rejected draw's target is overwritten by the
    /// next one). Then the swaps, in the same order as numpy. Same draws, same
    /// targets, same permutation, same stream position afterwards.
    pub fn shuffle<T: Copy>(&mut self, items: &mut [T]) {
        let n = items.len();
        if n < 2 {
            return;
        }
        assert!(n - 1 <= u32::MAX as usize, "shuffle bound exceeds 32 bits");
        let mut targets = std::mem::take(&mut self.targets);
        targets.clear();
        targets.resize(n, 0);
        let mut i = (n - 1) as u32;
        while i >= 1 {
            // Positions [lo, i] share the mask 2^(k+1) - 1.
            let lo = 1u32 << (31 - i.leading_zeros());
            let mask = lo.wrapping_shl(1).wrapping_sub(1) | lo;
            while i >= lo {
                if self.pos == self.len {
                    self.refill();
                }
                for &v in &self.buf[self.pos..self.len] {
                    let j = v & mask;
                    targets[i as usize] = j;
                    i -= (j <= i) as u32;
                    self.pos += 1;
                    if i < lo {
                        break;
                    }
                }
            }
        }
        for i in (1..n).rev() {
            items.swap(i, targets[i] as usize);
        }
        self.targets = targets;
    }

    /// Consume exactly the draws [`shuffle`](Self::shuffle) of `n` items
    /// would, without shuffling anything (only `random_interval`'s
    /// accept/reject process). Used to find where later columns' draws start.
    pub fn skip_shuffle(&mut self, n: usize) {
        if n < 2 {
            return;
        }
        assert!(n - 1 <= u32::MAX as usize, "shuffle bound exceeds 32 bits");
        let mut i = (n - 1) as u32;
        while i >= 1 {
            let lo = 1u32 << (31 - i.leading_zeros());
            let mask = lo.wrapping_shl(1).wrapping_sub(1) | lo;
            while i >= lo {
                if self.pos == self.len {
                    self.refill();
                }
                for &v in &self.buf[self.pos..self.len] {
                    i -= ((v & mask) <= i) as u32;
                    self.pos += 1;
                    if i < lo {
                        break;
                    }
                }
            }
        }
    }
}

/// Visit each column of a row-major `nrows x ncols` matrix after shuffling it
/// independently, in column order, with one generator seeded by `seed` —
/// exactly `np.apply_along_axis(default_rng(seed).permutation, 0, matrix)`.
///
/// `visit(j, column)` receives the permuted column `j` (length `nrows`).
pub fn for_each_permuted_column<T: Copy>(
    data: &[T],
    nrows: usize,
    ncols: usize,
    seed: u64,
    mut visit: impl FnMut(usize, &[T]),
) {
    let mut rng = U32Stream::from_seed(seed);
    let mut col: Vec<T> = Vec::with_capacity(nrows);
    for j in 0..ncols {
        col.clear();
        col.extend((0..nrows).map(|t| data[t * ncols + j]));
        rng.shuffle(&mut col);
        visit(j, &col);
    }
}

/// The permuted matrix itself (row-major), for callers that need it whole.
pub fn permute_columns<T: Copy + Default>(data: &[T], nrows: usize, ncols: usize, seed: u64) -> Vec<T> {
    let mut out = vec![T::default(); data.len()];
    for_each_permuted_column(data, nrows, ncols, seed, |j, col| {
        for (t, &v) in col.iter().enumerate() {
            out[t * ncols + j] = v;
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u32_stream_draws_match_next_u32() {
        for seed in [0u64, 1, 12345, 12345 + 100, u64::MAX] {
            let mut a = Pcg64::from_seed(seed);
            let mut b = U32Stream::from_seed(seed);
            for _ in 0..5000 {
                assert_eq!(a.next_u32(), b.next_u32());
            }
        }
        // A generator with a buffered half pending.
        let mut a = Pcg64::from_seed(7);
        a.next_u32();
        let mut b = U32Stream::new(a.clone());
        for _ in 0..3000 {
            assert_eq!(a.next_u32(), b.next_u32());
        }
    }

    #[test]
    fn advance_matches_stepping() {
        let mut a = Pcg64::from_seed(12345);
        let mut b = a.clone();
        for delta in [0u128, 1, 2, 3, 100, 1000, 12_345] {
            for _ in 0..delta {
                a.step();
            }
            b.advance(delta);
            assert_eq!(a.state, b.state, "delta {delta}");
        }
    }

    #[test]
    fn stream_at_position_continues_the_stream() {
        for seed in [12345u64, 12445, 7] {
            let mut full = U32Stream::from_seed(seed);
            let mut drawn = Vec::new();
            for _ in 0..3000 {
                drawn.push(full.next_u32());
            }
            for p in [0u64, 1, 2, 3, 511, 512, 513, 1023, 1024, 1025, 2999] {
                let mut s = U32Stream::at(seed, p);
                assert_eq!(s.position(), p);
                for (k, &want) in drawn[p as usize..].iter().enumerate().take(600) {
                    assert_eq!(s.next_u32(), want, "seed {seed}, position {p}, draw {k}");
                }
            }
        }
    }

    #[test]
    fn skip_shuffle_consumes_the_same_draws() {
        for seed in [12345u64, 999] {
            let mut shuffled = U32Stream::from_seed(seed);
            let mut skipped = U32Stream::from_seed(seed);
            for n in [0usize, 1, 2, 3, 255, 256, 257, 1000, 4097, 70_001, 5, 65_536] {
                let mut v: Vec<u32> = (0..n as u32).collect();
                shuffled.shuffle(&mut v);
                skipped.skip_shuffle(n);
                assert_eq!(skipped.position(), shuffled.position(), "seed {seed}, n {n}");
                // A stream started at that position continues identically.
                let mut resumed = U32Stream::at(seed, shuffled.position());
                let mut a: Vec<u32> = (0..300).collect();
                let mut b = a.clone();
                shuffled.clone_stream().shuffle(&mut a);
                resumed.shuffle(&mut b);
                assert_eq!(a, b, "seed {seed}, after n {n}");
            }
        }
    }

    #[test]
    fn u32_stream_shuffle_matches_reference() {
        let mut sizes: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 7, 8, 9, 100, 255, 256, 257, 511, 512, 513, 1000];
        for k in [10u32, 12, 16] {
            sizes.extend([(1usize << k) - 1, 1 << k, (1 << k) + 1]);
        }
        sizes.push(70_001);
        for seed in [12345u64, 12345 + 37, 999] {
            // One stream across many shuffles, like one generator across the
            // columns of an alignment: the carry-over between them must match.
            let mut reference = Pcg64::from_seed(seed);
            let mut stream = U32Stream::from_seed(seed);
            for (round, &n) in sizes.iter().chain(sizes.iter().rev()).enumerate() {
                let mut want: Vec<u32> = (0..n as u32).collect();
                let mut got = want.clone();
                reference.shuffle(&mut want);
                stream.shuffle(&mut got);
                assert_eq!(got, want, "seed {seed}, round {round}, n {n}");
            }
            assert_eq!(stream.next_u32(), reference.next_u32(), "stream position after shuffles");
        }
    }
}
