//! numpy-exact numerics (numpy 2.4.6).
//!
//! Floating-point results depend on the order of operations, so each helper
//! reproduces numpy's exact evaluation order:
//!
//! * [`pairwise_sum`]: reductions along a *contiguous* axis (`np.sum` of a 1-D
//!   array, `axis=-1` of a C-contiguous 2-D array, or the whole of a contiguous
//!   array) use numpy's pairwise summation (8 accumulators, block 128).
//! * reductions along a *non-contiguous* axis (e.g. `axis=0` of a C-contiguous
//!   (R, C) array with C > 1) accumulate row after row, i.e. sequentially;
//!   callers write those as plain loops. When C == 1 numpy drops the unit axis
//!   and the reduction becomes contiguous, so pairwise again.
//! * [`quantile_sorted`] ports `np.quantile(..., method="linear")` including the
//!   `_lerp` branch for gamma >= 0.5.

/// numpy `pairwise_sum` for float64 (`loops_utils.h.src`).
pub fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = -0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - (n % 8) {
            r[0] += a[i];
            r[1] += a[i + 1];
            r[2] += a[i + 2];
            r[3] += a[i + 3];
            r[4] += a[i + 4];
            r[5] += a[i + 5];
            r[6] += a[i + 6];
            r[7] += a[i + 7];
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// `np.sum` of a contiguous float64 array (as a reduction with identity 0).
pub fn sum(a: &[f64]) -> f64 {
    pairwise_sum(a)
}

/// `np.mean` of a contiguous float64 array.
pub fn mean(a: &[f64]) -> f64 {
    sum(a) / a.len() as f64
}

/// Sort ascending with NaN last (the order `np.partition`/`np.sort` imply).
pub fn sort_nan_last(v: &mut [f64]) {
    v.sort_unstable_by(|a, b| a.total_cmp(b));
    // total_cmp orders -NaN first; numpy puts every NaN last.
    let n_neg_nan = v.iter().take_while(|x| x.is_nan()).count();
    if n_neg_nan > 0 {
        v.rotate_left(n_neg_nan);
    }
}

/// `np.median` of a 1-D array (consumes a scratch copy).
pub fn median(a: &[f64]) -> f64 {
    let mut v = a.to_vec();
    median_in_place(&mut v)
}

/// `np.median` using `v` as scratch space (it is reordered).
pub fn median_in_place(v: &mut [f64]) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    // Any NaN makes the median NaN (`_median_nancheck`).
    if v.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }
    let half = n / 2;
    if n % 2 == 1 {
        let (_, m, _) = v.select_nth_unstable_by(half, |a, b| a.total_cmp(b));
        // mean of one element: x / 1
        (-0.0 + *m) / 1.0
    } else {
        let (lower, hi, _) = v.select_nth_unstable_by(half, |a, b| a.total_cmp(b));
        let hi = *hi;
        let lo = lower.iter().copied().fold(f64::NEG_INFINITY, |m, x| if x.total_cmp(&m).is_gt() { x } else { m });
        ((-0.0 + lo) + hi) / 2.0
    }
}

/// k-th smallest value (`np.partition(a, k)[k]`), NaN sorting last.
pub fn kth_smallest(a: &[f64], k: usize) -> f64 {
    let mut v = a.to_vec();
    sort_nan_last(&mut v);
    v[k]
}

/// `np.linspace(start, stop, num)` for float64 scalars (endpoint=True).
pub fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    if num == 0 {
        return Vec::new();
    }
    let div = num - 1;
    let delta = stop - start;
    let mut y: Vec<f64> = (0..num).map(|i| i as f64).collect();
    if div > 0 {
        let step = delta / div as f64;
        if step == 0.0 {
            for v in y.iter_mut() {
                *v = *v / div as f64 * delta;
            }
        } else {
            for v in y.iter_mut() {
                *v *= step;
            }
        }
    } else {
        for v in y.iter_mut() {
            *v *= delta;
        }
    }
    for v in y.iter_mut() {
        *v += start;
    }
    if num > 1 {
        y[num - 1] = stop;
    }
    y
}

/// `np.linspace(0, 1, K + 2)[1:-1]`: the interior quantile positions.
pub fn interior_positions(k: usize) -> Vec<f64> {
    let full = linspace(0.0, 1.0, k + 2);
    full[1..k + 1].to_vec()
}

/// `np.linspace(0, L, n + 1, dtype=int)` (chunk boundaries).
pub fn linspace_int(stop: usize, num: usize) -> Vec<usize> {
    linspace(0.0, stop as f64, num).into_iter().map(|v| v.floor() as usize).collect()
}

/// Precomputed index/gamma pairs for `np.quantile(x, q, method="linear")` on
/// arrays of a fixed length `n`.
#[derive(Clone, Debug)]
pub struct QuantilePlan {
    pub n: usize,
    prev: Vec<usize>,
    next: Vec<usize>,
    gamma: Vec<f64>,
    /// The distinct sorted positions the plan reads, ascending.
    ranks: Vec<usize>,
    /// Index into `ranks` of `prev[i]` / `next[i]`.
    prev_pos: Vec<usize>,
    next_pos: Vec<usize>,
}

impl QuantilePlan {
    pub fn new(n: usize, quantiles: &[f64]) -> Self {
        let mut prev = Vec::with_capacity(quantiles.len());
        let mut next = Vec::with_capacity(quantiles.len());
        let mut gamma = Vec::with_capacity(quantiles.len());
        for &q in quantiles {
            let virtual_index = (n as f64 - 1.0) * q;
            let mut p = virtual_index.floor();
            let mut nx = p + 1.0;
            if virtual_index >= n as f64 - 1.0 {
                p = -1.0;
                nx = -1.0;
            }
            if virtual_index < 0.0 {
                p = 0.0;
                nx = 0.0;
            }
            if virtual_index.is_nan() {
                p = -1.0;
                nx = -1.0;
            }
            let pi = p as i64;
            let ni = nx as i64;
            // gamma uses the integer previous index (after the bounds fixes)
            gamma.push(virtual_index - pi as f64);
            let wrap = |i: i64| if i < 0 { (n as i64 + i) as usize } else { i as usize };
            prev.push(wrap(pi));
            next.push(wrap(ni));
        }
        let mut ranks: Vec<usize> = prev.iter().chain(&next).copied().collect();
        ranks.sort_unstable();
        ranks.dedup();
        let pos = |r: &usize| ranks.binary_search(r).expect("rank listed");
        let prev_pos = prev.iter().map(pos).collect();
        let next_pos = next.iter().map(pos).collect();
        QuantilePlan { n, prev, next, gamma, ranks, prev_pos, next_pos }
    }

    /// The distinct sorted positions the quantiles read, ascending.
    pub fn ranks(&self) -> &[usize] {
        &self.ranks
    }

    /// Quantiles from the values at [`ranks`](Self::ranks) of the sorted array
    /// (`at_ranks[i]` is the value at `ranks()[i]`; no NaN in the array).
    /// Same arithmetic as [`apply_sorted`](Self::apply_sorted).
    pub fn apply_ranked(&self, at_ranks: &[f64], out: &mut [f64]) {
        debug_assert_eq!(at_ranks.len(), self.ranks.len());
        for i in 0..self.gamma.len() {
            out[i] = lerp(at_ranks[self.prev_pos[i]], at_ranks[self.next_pos[i]], self.gamma[i]);
        }
    }

    pub fn len(&self) -> usize {
        self.gamma.len()
    }

    pub fn is_empty(&self) -> bool {
        self.gamma.is_empty()
    }

    /// Quantiles of `sorted` (ascending, NaN last), written into `out`.
    pub fn apply_sorted(&self, sorted: &[f64], out: &mut [f64]) {
        debug_assert_eq!(sorted.len(), self.n);
        if sorted.last().is_some_and(|x| x.is_nan()) {
            out.fill(f64::NAN);
            return;
        }
        for i in 0..self.gamma.len() {
            out[i] = lerp(sorted[self.prev[i]], sorted[self.next[i]], self.gamma[i]);
        }
    }
}

/// numpy `_lerp(a, b, t)`.
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    let diff = b - a;
    if t >= 0.5 { b - diff * (1.0 - t) } else { a + diff * t }
}

/// `np.quantile(a, quantiles)` (method="linear") for a 1-D array.
pub fn quantile(a: &[f64], quantiles: &[f64]) -> Vec<f64> {
    let mut v = a.to_vec();
    sort_nan_last(&mut v);
    let plan = QuantilePlan::new(v.len(), quantiles);
    let mut out = vec![0.0; quantiles.len()];
    plan.apply_sorted(&v, &mut out);
    out
}

/// `np.percentile(a, p)` for a scalar percentile.
pub fn percentile(a: &[f64], p: f64) -> f64 {
    quantile(a, &[p / 100.0])[0]
}

/// Below this length (or when targets are dense) [`select_ranks`] just sorts.
const SELECT_SMALL: usize = 64;
/// Sort instead of bucketing when `ranks * SELECT_DENSE >= len`.
const SELECT_DENSE: usize = 8;
/// About `len / SELECT_BUCKET_DIV` buckets per level (power of two, 16..4096).
const SELECT_BUCKET_DIV: usize = 4;
/// Bucketing levels before falling back to a sort (each level shrinks the key
/// range by at least 2^3, so this is rarely reached).
const SELECT_MAX_DEPTH: usize = 12;

/// Reusable buffers for [`select_ranks`].
#[derive(Default)]
pub struct SelectScratch {
    keys: Vec<u64>,
    out: Vec<u64>,
    levels: Vec<Vec<u64>>,
    counts: Vec<u32>,
}

/// Order-preserving map from `f64` to `u64`: `a.total_cmp(b) == key(a).cmp(key(b))`
/// (negative values have all bits flipped, non-negative ones the sign bit).
#[inline]
pub fn total_order_key(x: f64) -> u64 {
    let b = x.to_bits();
    if b >> 63 == 1 { !b } else { b | (1 << 63) }
}

/// Inverse of [`total_order_key`].
#[inline]
pub fn from_total_order_key(k: u64) -> f64 {
    f64::from_bits(if k >> 63 == 1 { k & !(1 << 63) } else { !k })
}

/// The values at sorted positions `ranks` (ascending, distinct, each <
/// `v.len()`) of `v` in `f64::total_cmp` order, i.e. `sorted(v)[r]` for each r,
/// written to `out`. `v` must not contain NaN.
///
/// Exact: values are mapped to their total-order keys, bucketed by the high
/// bits of `key - min_key` (a monotone map), only buckets that contain a
/// wanted position are kept, and those are refined the same way or sorted.
/// Equal keys are identical bit patterns, and the value at a position of a
/// total order does not depend on how it is found, so the result is
/// bit-for-bit what a full sort gives. Bucketing on the key rather than on
/// the value spreads skewed data well (the exponent is in the high bits).
pub fn select_ranks(v: &[f64], ranks: &[usize], out: &mut [f64], scratch: &mut SelectScratch) {
    select_ranks_with(v, ranks, out, scratch, SELECT_SMALL, SELECT_DENSE);
}

/// [`select_ranks`] with explicit thresholds (tests force the bucketing path).
pub fn select_ranks_with(
    v: &[f64],
    ranks: &[usize],
    out: &mut [f64],
    scratch: &mut SelectScratch,
    small: usize,
    dense: usize,
) {
    debug_assert!(ranks.windows(2).all(|w| w[0] < w[1]));
    debug_assert!(ranks.last().is_none_or(|&r| r < v.len()));
    debug_assert!(!v.iter().any(|x| x.is_nan()));
    let SelectScratch { keys, out: kout, levels, counts } = scratch;
    keys.clear();
    keys.extend(v.iter().map(|&x| total_order_key(x)));
    kout.clear();
    kout.resize(ranks.len(), 0);
    select_keys(keys, ranks, kout, levels, counts, 0, small, dense);
    for (o, &k) in out.iter_mut().zip(kout.iter()) {
        *o = from_total_order_key(k);
    }
}

/// Keys at `ranks` of a small slice: one rank, or two adjacent ones (the
/// usual `prev`/`next` pair of a quantile), by partial selection; otherwise a
/// sort. Keys are totally ordered integers, so all three are exact.
fn sort_and_read(v: &mut [u64], ranks: &[usize], out: &mut [u64]) {
    match *ranks {
        [r] => out[0] = *v.select_nth_unstable(r).1,
        [r, r1] if r1 == r + 1 => {
            let (_, m, right) = v.select_nth_unstable(r);
            out[0] = *m;
            out[1] = *right.iter().min().expect("r + 1 < len");
        }
        _ => {
            v.sort_unstable();
            for (o, &r) in out.iter_mut().zip(ranks) {
                *o = v[r];
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn select_keys(
    v: &mut [u64],
    ranks: &[usize],
    out: &mut [u64],
    levels: &mut Vec<Vec<u64>>,
    counts: &mut Vec<u32>,
    depth: usize,
    small: usize,
    dense: usize,
) {
    if ranks.is_empty() {
        return;
    }
    let n = v.len();
    if n <= small || ranks.len() * dense >= n || depth >= SELECT_MAX_DEPTH {
        return sort_and_read(v, ranks, out);
    }
    let (mut lo, mut hi) = (v[0], v[0]);
    for &k in v.iter() {
        lo = lo.min(k);
        hi = hi.max(k);
    }
    if lo == hi {
        out.fill(lo);
        return;
    }
    // Bucket = (key - lo) >> shift, with at most 2^nb_bits buckets. lo lands
    // in bucket 0 and hi in the last one (>= 1), so every bucket is smaller
    // than v and each level shrinks the key range by >= 2^(nb_bits - 1).
    let range = hi - lo;
    let nb_bits = ((n / SELECT_BUCKET_DIV).max(16).ilog2()).min(12);
    let shift = (64 - range.leading_zeros()).saturating_sub(nb_bits);
    let nb = (range >> shift) as usize + 1;
    counts.clear();
    counts.resize(nb, 0);
    for &k in v.iter() {
        counts[((k - lo) >> shift) as usize] += 1;
    }
    // Buckets holding a wanted position: (bucket, first position, size, ranks).
    let mut needed: Vec<(usize, usize, usize, std::ops::Range<usize>)> = Vec::new();
    let (mut start, mut ri) = (0usize, 0usize);
    for (b, &c) in counts.iter().enumerate() {
        let end = start + c as usize;
        let r0 = ri;
        while ri < ranks.len() && ranks[ri] < end {
            ri += 1;
        }
        if ri > r0 {
            needed.push((b, start, c as usize, r0..ri));
        }
        start = end;
        if ri == ranks.len() {
            break;
        }
    }
    // Branch-free gather of the needed buckets: each bucket gets a write
    // cursor and a step; buckets that are not needed all write to one dummy
    // slot at the end (step 0), so no data-dependent branch is taken.
    let total: usize = needed.iter().map(|nd| nd.2).sum();
    let dummy = total as u32;
    counts.fill(dummy);
    let mut steps = vec![0u32; nb];
    let mut cursor = 0u32;
    for nd in &needed {
        counts[nd.0] = cursor;
        steps[nd.0] = 1;
        cursor += nd.2 as u32;
    }
    if levels.len() <= depth {
        levels.resize_with(depth + 1, Vec::new);
    }
    let mut buf = std::mem::take(&mut levels[depth]);
    buf.clear();
    buf.resize(total + 1, 0);
    for &k in v.iter() {
        let b = ((k - lo) >> shift) as usize;
        let c = counts[b];
        buf[c as usize] = k;
        counts[b] = c + steps[b];
    }
    buf.truncate(total);
    let mut off = 0usize;
    let mut rel: Vec<usize> = Vec::new();
    for (_, first, size, rr) in &needed {
        rel.clear();
        rel.extend(ranks[rr.clone()].iter().map(|r| r - first));
        select_keys(&mut buf[off..off + size], &rel, &mut out[rr.clone()], levels, counts, depth + 1, small, dense);
        off += size;
    }
    levels[depth] = buf;
}

#[cfg(test)]
mod tests {
    use super::*;

    struct XorShift(u64);
    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Value sets that stress the bucketing: spreads, ties, signed zeros,
    /// infinities, subnormals, constants and heavy tails.
    fn cases() -> Vec<Vec<f64>> {
        let mut r = XorShift(0x853c_49e6_748f_ea9b);
        let mut out = Vec::new();
        out.push((0..5000).map(|_| r.unit() * 100.0 - 20.0).collect());
        out.push((0..5000).map(|_| (r.next() % 7) as f64).collect());
        out.push((0..3000).map(|i| if i % 2 == 0 { -0.0 } else { 0.0 }).collect());
        out.push(
            (0..3000)
                .map(|i| match i % 5 {
                    0 => -0.0,
                    1 => 0.0,
                    2 => r.unit() - 0.5,
                    3 => 5e-324 * (r.next() % 4) as f64,
                    _ => -5e-324,
                })
                .collect(),
        );
        out.push(
            (0..2000)
                .map(|i| match i % 97 {
                    0 => f64::INFINITY,
                    1 => f64::NEG_INFINITY,
                    _ => r.unit(),
                })
                .collect(),
        );
        out.push(vec![3.25; 1000]);
        // Exponential tail: each bucketing level only peels off the top.
        out.push((0..1500).map(|i| 2f64.powi(i % 900)).collect());
        // Bulk plus far outliers (like biased taxa).
        out.push(
            (0..8000)
                .map(|i| if i % 1000 == 7 { 1e9 * r.unit() } else { 10.0 + r.unit() })
                .collect(),
        );
        out.push((0..4000).map(|_| f64::MAX * (r.unit() - 0.5)).collect());
        out.push(vec![1.0]);
        out
    }

    #[test]
    fn select_ranks_matches_full_sort() {
        let mut r = XorShift(0x9e37_79b9_7f4a_7c15);
        let mut scratch = SelectScratch::default();
        for (ci, v) in cases().into_iter().enumerate() {
            let n = v.len();
            let mut sorted = v.clone();
            sorted.sort_unstable_by(f64::total_cmp);
            let mut rank_sets: Vec<Vec<usize>> = vec![
                (0..n).collect(),
                vec![0],
                vec![n - 1],
                vec![0, n - 1],
            ];
            for _ in 0..4 {
                let mut rs: Vec<usize> = (0..1 + r.next() as usize % 400).map(|_| r.next() as usize % n).collect();
                rs.sort_unstable();
                rs.dedup();
                rank_sets.push(rs);
            }
            rank_sets.push(QuantilePlan::new(n, &interior_positions(200.min(n))).ranks().to_vec());
            for rs in rank_sets.iter_mut() {
                rs.dedup(); // e.g. [0, n - 1] when n == 1
            }
            for ranks in &rank_sets {
                for &(small, dense) in &[(2usize, 0usize), (16, 1), (SELECT_SMALL, SELECT_DENSE)] {
                    let mut got = vec![f64::NAN; ranks.len()];
                    select_ranks_with(&v, ranks, &mut got, &mut scratch, small, dense);
                    for (g, &rk) in got.iter().zip(ranks) {
                        assert_eq!(
                            g.to_bits(),
                            sorted[rk].to_bits(),
                            "case {ci}, rank {rk}, small {small}, dense {dense}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn total_order_key_matches_total_cmp() {
        let vals = [
            f64::NEG_INFINITY, -f64::MAX, -1.5, -f64::MIN_POSITIVE, -5e-324, -0.0, 0.0, 5e-324,
            f64::MIN_POSITIVE, 1.0, 1.5, f64::MAX, f64::INFINITY, f64::NAN, -f64::NAN,
        ];
        for &a in &vals {
            assert_eq!(from_total_order_key(total_order_key(a)).to_bits(), a.to_bits());
            for &b in &vals {
                assert_eq!(total_order_key(a).cmp(&total_order_key(b)), a.total_cmp(&b), "{a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn apply_ranked_matches_apply_sorted() {
        let mut r = XorShift(42);
        for &n in &[1usize, 2, 3, 50, 199, 200, 201, 1000, 12944] {
            let mut sorted: Vec<f64> = (0..n).map(|_| r.unit() * 10.0).collect();
            sorted.sort_unstable_by(f64::total_cmp);
            for k in [1usize, 7, 200] {
                let plan = QuantilePlan::new(n, &interior_positions(k));
                let mut want = vec![0.0; k];
                plan.apply_sorted(&sorted, &mut want);
                let at: Vec<f64> = plan.ranks().iter().map(|&i| sorted[i]).collect();
                let mut got = vec![0.0; k];
                plan.apply_ranked(&at, &mut got);
                for (g, w) in got.iter().zip(&want) {
                    assert_eq!(g.to_bits(), w.to_bits(), "n {n}, k {k}");
                }
            }
        }
    }
}
