//! Permutation null distribution (`PermutationTest` in Python) and the
//! derived null state used by the pruner.

use rayon::prelude::*;

use crate::chi2::{self, Encoded};
use crate::npy;
use crate::rng::U32Stream;

/// Base seed used by WitChi: permutation `i` uses `12345 + i`.
pub const SEED_BASE: u64 = 12345;

/// Visit the columns of a column-major encoded alignment after shuffling each
/// one with a single generator seeded by `seed`, in column order
/// (`np.apply_along_axis(default_rng(seed).permutation, 0, a)`).
pub fn for_each_permuted_column(enc: &Encoded, seed: u64, mut visit: impl FnMut(usize, &[i8])) {
    let mut rng = U32Stream::from_seed(seed);
    let mut col = vec![0i8; enc.nrows];
    for j in 0..enc.ncols {
        col.copy_from_slice(enc.column(j));
        rng.shuffle(&mut col);
        visit(j, &col);
    }
}

/// The whole permuted alignment (column-major). The delta null streams the
/// same columns instead; this is the reference it is tested against.
pub fn permuted(enc: &Encoded, seed: u64) -> Encoded {
    let mut codes = vec![0i8; enc.nrows * enc.ncols];
    let t = enc.nrows;
    for_each_permuted_column(enc, seed, |j, col| codes[j * t..(j + 1) * t].copy_from_slice(col));
    Encoded::from_codes(codes, enc.nrows, enc.ncols)
}

/// Columns shuffled and counted together by [`permuted_counts`].
pub(crate) const COUNT_BATCH: usize = 32;

/// Column segments of one permutation that can be processed in parallel:
/// segment k covers columns `starts[k]..starts[k + 1]` (the last runs to the
/// end), and its shuffles start at 32-bit draw `positions[k]` of the
/// permutation's stream.
///
/// Each column's shuffle consumes a number of draws that depends only on the
/// draws themselves and the column length, never on the residues, so the
/// positions are found by a cheap draws-only pass
/// ([`U32Stream::skip_shuffle`]) and each segment then jumps its own stream
/// there ([`U32Stream::at`]). Every column gets exactly the draws it gets in
/// one sequential pass, so the permutation is the same.
#[derive(Clone, Debug)]
pub struct Segments {
    pub starts: Vec<usize>,
    pub positions: Vec<u64>,
}

impl Segments {
    /// One segment: the whole alignment from draw 0 (no pre-pass).
    pub fn whole() -> Self {
        Segments { starts: vec![0], positions: vec![0] }
    }

    /// Segments starting at the given column indices (ascending, first 0).
    pub fn at_columns(enc: &Encoded, seed: u64, starts: Vec<usize>) -> Self {
        debug_assert_eq!(starts.first(), Some(&0));
        let mut positions = vec![0u64];
        if starts.len() > 1 {
            let mut stream = U32Stream::from_seed(seed);
            for k in 1..starts.len() {
                for _ in starts[k - 1]..starts[k] {
                    stream.skip_shuffle(enc.nrows);
                }
                positions.push(stream.position());
            }
        }
        Segments { starts, positions }
    }

    pub fn len(&self) -> usize {
        self.starts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    /// Columns of segment `k`.
    pub fn range(&self, k: usize, ncols: usize) -> std::ops::Range<usize> {
        self.starts[k]..self.starts.get(k + 1).copied().unwrap_or(ncols)
    }
}

/// Segments per permutation so that `permutations` permutations keep the
/// current thread pool busy: the count minimizing the estimated makespan
/// `ceil(P * n / threads) * cost / n + skip` (skip only when n > 1), where
/// `cost` and `skip` are the per-cell costs of the stage and of the draws-only
/// pass. With enough permutations per thread this is 1 (no pre-pass).
pub fn segment_count(permutations: usize, units: usize, cost: f64, skip: f64) -> usize {
    let threads = rayon::current_num_threads().max(1);
    let mut best = (f64::INFINITY, 1);
    for n in 1..=units.clamp(1, 64) {
        let rounds = (permutations * n).div_ceil(threads) as f64;
        let est = rounds * cost / n as f64 + if n > 1 { skip } else { 0.0 };
        if est < best.0 * (1.0 - 1e-9) {
            best = (est, n);
        }
    }
    best.1
}

/// Relative per-cell costs (ns) used by [`segment_count`]: the draws-only
/// pass, and shuffling plus counting.
pub(crate) const SKIP_COST: f64 = 2.5;
pub(crate) const COUNT_COST: f64 = 6.0;

/// Column starts of `nseg` segments made of whole `unit`-column units.
pub fn segment_starts(ncols: usize, unit: usize, nseg: usize) -> Vec<usize> {
    let units = ncols.div_ceil(unit).max(1);
    let nseg = nseg.clamp(1, units);
    (0..nseg).map(|k| (k * units / nseg) * unit).collect()
}

/// Raw (C, T) residue counts of the alignment permuted with `seed` (the
/// columns of `np.apply_along_axis(default_rng(seed).permutation, 0, a)`),
/// without materializing it.
pub fn permuted_counts(enc: &Encoded, n_chars: usize, seed: u64) -> Vec<i64> {
    permuted_counts_in(enc, n_chars, seed, &Segments::whole())
}

/// [`permuted_counts`] with the columns split into `segs`, counted in
/// parallel.
///
/// Within a segment, columns are shuffled in order with the segment's
/// stream, a batch at a time, and each batch is added to a taxon-major
/// (T, C + 1) table in one sequential sweep (slot 0 takes all negative codes,
/// i.e. residues outside the char set). The segments' tables are summed.
/// Counts are integers, so the result equals counting column by column.
pub fn permuted_counts_in(enc: &Encoded, n_chars: usize, seed: u64, segs: &Segments) -> Vec<i64> {
    let t_n = enc.nrows;
    let w = n_chars + 1;
    let mut counts = vec![0i64; n_chars * t_n];
    if t_n == 0 {
        return counts;
    }
    let count_segment = |k: usize| count_columns(enc, n_chars, seed, segs.positions[k], segs.range(k, enc.ncols));
    let table = if segs.len() == 1 {
        count_segment(0)
    } else {
        (0..segs.len())
            .into_par_iter()
            .map(count_segment)
            .reduce_with(|mut a, b| {
                for (x, y) in a.iter_mut().zip(&b) {
                    *x += y;
                }
                a
            })
            .expect("at least one segment")
    };
    for t in 0..t_n {
        for c in 0..n_chars {
            counts[c * t_n + t] = table[t * w + c + 1] as i64;
        }
    }
    counts
}

/// Taxon-major (T, C + 1) residue counts of columns `cols` shuffled with the
/// stream of `seed` from draw `position`. A plain function rather than a
/// closure: captured-by-reference sizes in the hot loop cost about 40 %.
fn count_columns(enc: &Encoded, n_chars: usize, seed: u64, position: u64, cols: std::ops::Range<usize>) -> Vec<u32> {
    let t_n = enc.nrows;
    let w = n_chars + 1;
    let mut rng = U32Stream::at(seed, position);
    let mut batch = vec![0i8; COUNT_BATCH.min(cols.len()) * t_n];
    let mut table = vec![0u32; t_n * w];
    let mut j = cols.start;
    while j < cols.end {
        let nb = COUNT_BATCH.min(cols.end - j);
        let batch = &mut batch[..nb * t_n];
        for (k, slot) in batch.chunks_exact_mut(t_n).enumerate() {
            slot.copy_from_slice(enc.column(j + k));
            rng.shuffle(slot);
        }
        for t in 0..t_n {
            let row = &mut table[t * w..(t + 1) * w];
            for k in 0..nb {
                row[(batch[k * t_n + t].max(-1) as isize + 1) as usize] += 1;
            }
        }
        j += nb;
    }
    table
}

/// Per-taxon chi^2 of one permuted alignment.
pub fn permuted_row_chi2(enc: &Encoded, n_chars: usize, seed: u64) -> Vec<f64> {
    permuted_row_chi2_in(enc, n_chars, seed, &Segments::whole())
}

fn permuted_row_chi2_in(enc: &Encoded, n_chars: usize, seed: u64, segs: &Segments) -> Vec<f64> {
    let counts = chi2::fudge_counts(permuted_counts_in(enc, n_chars, seed, segs));
    let e = chi2::expected_observed(&counts, n_chars, enc.nrows);
    chi2::row_chi2(&e, &counts, n_chars, enc.nrows)
}

/// The (P, T) matrix of permuted per-taxon chi^2 and its summaries
/// (`PermutationTest.compute_null` / `_summarize_null`).
#[derive(Clone, Debug)]
pub struct NullMatrix {
    pub permutations: usize,
    pub ntaxa: usize,
    /// Row-major (P, T); also the pooled null (`ravel()`).
    pub pooled: Vec<f64>,
    /// Per-permutation alignment chi^2 (`np.sum(matrix, axis=1)`).
    pub sums: Vec<f64>,
}

impl NullMatrix {
    /// Runs on the current rayon pool (results do not depend on its size).
    pub fn compute(enc: &Encoded, n_chars: usize, permutations: usize) -> Self {
        // With fewer permutations than the pool can keep busy, each one is
        // also split into column segments (same counts, see `Segments`).
        let nseg = segment_count(permutations, enc.ncols.div_ceil(COUNT_BATCH), COUNT_COST, SKIP_COST);
        let starts = segment_starts(enc.ncols, COUNT_BATCH, nseg);
        let rows: Vec<Vec<f64>> = (0..permutations)
            .into_par_iter()
            .map(|i| {
                let seed = SEED_BASE + i as u64;
                let segs = Segments::at_columns(enc, seed, starts.clone());
                permuted_row_chi2_in(enc, n_chars, seed, &segs)
            })
            .collect();
        let sums = rows.iter().map(|r| npy::sum(r)).collect();
        let pooled = rows.concat();
        NullMatrix { permutations, ntaxa: enc.nrows, pooled, sums }
    }
}

/// Run `f` on a rayon pool with `threads` workers (1 = current thread order
/// is irrelevant: results are collected in index order either way).
pub fn with_threads<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("thread pool")
        .install(f)
}

/// Alignment-level empirical p-value: `#(null >= x) / len(null)`.
pub fn empirical_p(x: f64, null: &[f64]) -> f64 {
    let n = null.iter().filter(|&&v| v >= x).count();
    n as f64 / null.len() as f64
}

/// A null distribution prepared for many `#(null >= x)` queries: the non-NaN
/// values sorted ascending (NaN never satisfies `>=`), plus the full length.
pub struct SortedNull {
    sorted: Vec<f64>,
    len: usize,
}

impl SortedNull {
    pub fn new(null: &[f64]) -> Self {
        let mut sorted: Vec<f64> = null.iter().copied().filter(|v| !v.is_nan()).collect();
        sorted.par_sort_unstable_by(|a, b| a.total_cmp(b));
        SortedNull { sorted, len: null.len() }
    }

    /// `#(null >= x)`, exactly as the elementwise comparison counts it.
    ///
    /// `v < x` is monotone over the sorted values (total order differs from IEEE
    /// only in placing -0.0 before 0.0, which compare equal either way).
    pub fn count_ge(&self, x: f64) -> usize {
        if x.is_nan() {
            return 0;
        }
        self.sorted.len() - self.sorted.partition_point(|&v| v < x)
    }

    /// Per-taxon Bonferroni empirical p-values, clipped to [0, 1].
    pub fn p_per_taxon(&self, per_row: &[f64]) -> Vec<f64> {
        let n_taxa = per_row.len() as f64;
        per_row
            .par_iter()
            .map(|&x| {
                let p = (self.count_ge(x) as f64 / self.len as f64) * n_taxa;
                p.max(0.0).min(1.0)
            })
            .collect()
    }
}

/// Per-taxon Bonferroni empirical p-values, clipped to [0, 1].
pub fn empirical_p_per_taxon(per_row: &[f64], null: &[f64]) -> Vec<f64> {
    SortedNull::new(null).p_per_taxon(per_row)
}

/// Robust Z-score: mean centring, MAD scaling (`utils._robust_zscore`).
pub fn robust_zscore(observed: f64, null_pool: &[f64]) -> f64 {
    let stats = RobustStats::of(null_pool);
    stats.zscore(observed)
}

/// Mean / median-absolute-deviation of a null pool.
#[derive(Clone, Copy, Debug)]
pub struct RobustStats {
    pub mean: f64,
    pub mad: f64,
}

impl RobustStats {
    pub fn of(null_pool: &[f64]) -> Self {
        let mean = npy::mean(null_pool);
        let med = npy::median(null_pool);
        let mut dev: Vec<f64> = null_pool.iter().map(|&v| (v - med).abs()).collect();
        let mad = npy::median_in_place(&mut dev);
        RobustStats { mean, mad }
    }

    pub fn zscore(&self, observed: f64) -> f64 {
        if self.mad == 0.0 {
            if observed == self.mean {
                0.0
            } else {
                let d = observed - self.mean;
                // np.sign(d) * inf (np.sign(nan) = nan, np.sign(0) = 0 -> nan)
                if d.is_nan() || d == 0.0 { f64::NAN } else { d.signum() * f64::INFINITY }
            }
        } else {
            (observed - self.mean) / (self.mad / 0.6745)
        }
    }
}

/// Null-derived quantities the pruner needs (`_compute_null_state`).
#[derive(Clone, Debug)]
pub struct NullState {
    pub null: NullMatrix,
    pub z_mean: f64,
    pub z_scale: f64,
    pub z_quantiles: Vec<f64>,
    pub biased_threshold: f64,
}

impl NullState {
    pub fn new(null: NullMatrix) -> Self {
        let pooled = &null.pooled;
        let stats = RobustStats::of(pooled);
        let z_mean = stats.mean;
        let z_scale = stats.mad / 0.6745;
        let null_z: Vec<f64> = pooled.iter().map(|&v| (v - z_mean) / z_scale).collect();
        let k = null_z.len().min(200);
        let z_quantiles = npy::quantile(&null_z, &npy::interior_positions(k));
        let biased_threshold = biased_taxa_threshold(pooled, null.ntaxa);
        NullState { null, z_mean, z_scale, z_quantiles, biased_threshold }
    }
}

/// `_compute_biased_taxa_threshold`: the (M - k - 1)-th order statistic of the
/// pooled null with k = floor(0.05 * M / n_taxa).
pub fn biased_taxa_threshold(pooled: &[f64], n_taxa: usize) -> f64 {
    let m = pooled.len();
    if m == 0 || n_taxa == 0 {
        return f64::INFINITY;
    }
    let k = (0.05 * m as f64 / n_taxa as f64).floor() as i64;
    let kth = m as i64 - k - 1;
    if kth < 0 {
        return f64::INFINITY;
    }
    npy::kth_smallest(pooled, kth as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permuted_counts_match_column_by_column_reference() {
        use crate::rng::Pcg64;
        let n_chars = 20;
        for &(t_n, ncols) in &[(1usize, 5usize), (7, 31), (50, 32), (50, 33), (333, 100), (1000, 65)] {
            let mut state = 0x853c_49e6_748f_ea9bu64 ^ (t_n * 31 + ncols) as u64;
            let codes: Vec<i8> = (0..t_n * ncols)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    let r = (state >> 40) % 23;
                    if r >= 20 { -1 } else { r as i8 }
                })
                .collect();
            let enc = Encoded::from_codes(codes, t_n, ncols);
            for seed in [SEED_BASE, SEED_BASE + 7, SEED_BASE + 1000] {
                // Reference: numpy's draw-by-draw shuffle, one column at a time.
                let mut rng = Pcg64::from_seed(seed);
                let mut want = vec![0i64; n_chars * t_n];
                let mut col = vec![0i8; t_n];
                for j in 0..ncols {
                    col.copy_from_slice(enc.column(j));
                    rng.shuffle(&mut col);
                    chi2::add_column_counts(&mut want, &col, t_n);
                }
                assert_eq!(permuted_counts(&enc, n_chars, seed), want, "T {t_n}, L {ncols}, seed {seed}");
            }
        }
    }

    #[test]
    fn segmented_counts_match_whole() {
        let n_chars = 20;
        let (t_n, ncols) = (257usize, 211usize);
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let codes: Vec<i8> = (0..t_n * ncols)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let r = (state >> 40) % 23;
                if r >= 20 { -1 } else { r as i8 }
            })
            .collect();
        let enc = Encoded::from_codes(codes, t_n, ncols);
        for seed in [SEED_BASE, SEED_BASE + 3] {
            let whole = permuted_counts(&enc, n_chars, seed);
            for starts in [vec![0], vec![0, 1], vec![0, 32, 64, 200], (0..ncols).collect::<Vec<_>>()] {
                let segs = Segments::at_columns(&enc, seed, starts.clone());
                assert_eq!(permuted_counts_in(&enc, n_chars, seed, &segs), whole, "starts {starts:?}");
            }
        }
        assert_eq!(segment_starts(100, 32, 3), vec![0, 32, 64]);
        assert_eq!(segment_starts(100, 32, 9), vec![0, 32, 64, 96]);
        // Plenty of permutations per thread: no segmentation.
        assert_eq!(segment_count(64 * rayon::current_num_threads(), 100, 50.0, 2.5), 1);
    }

    #[test]
    fn sorted_null_counts_like_elementwise_ge() {
        let mut null = vec![
            f64::NAN, -0.0, 0.0, 1.0, 1.0, 1.0, -2.5, f64::INFINITY, f64::NEG_INFINITY,
            -f64::NAN, 3.25, 1e-300, -1e-300, 7.0,
        ];
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for _ in 0..500 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            null.push(((state % 41) as f64 - 20.0) / 4.0);
        }
        let sn = SortedNull::new(&null);
        let mut probes = null.clone();
        probes.extend([0.5, -100.0, 100.0, f64::MAX, f64::MIN, f64::NAN]);
        for &x in &probes {
            let want = null.iter().filter(|&&v| v >= x).count();
            assert_eq!(sn.count_ge(x), want, "x = {x:?}");
        }
        let per_row = [1.0, -0.0, f64::NAN, 2.0];
        let want: Vec<f64> = per_row
            .iter()
            .map(|&x| {
                let p = (null.iter().filter(|&&v| v >= x).count() as f64 / null.len() as f64) * 4.0;
                p.max(0.0).min(1.0)
            })
            .collect();
        let got = sn.p_per_taxon(&per_row);
        for (g, w) in got.iter().zip(&want) {
            assert_eq!(g.to_bits(), w.to_bits());
        }
    }
}
