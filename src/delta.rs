//! Per-column delta scores: how much the active score drops when column j is
//! removed (`ChiSquareCalculator._reduce_chi2_per_col_removed` and the
//! `calculate_*_difference` methods), plus the delta-null distribution
//! (`AlignmentPruner._compute_null_deltas`).
//!
//! The column range is cut into the same chunks and L2-sized blocks as the
//! Python code, because block width changes which BLAS routine runs and how
//! numpy orders the axis-0 reductions (see `npy` module docs). Blocks are then
//! processed in parallel; each block's result depends only on its own columns.

use rayon::prelude::*;

use crate::blas;
use crate::chi2::{self, Encoded};
use crate::npy::{self, QuantilePlan};
use crate::permutation::{self, NullState, SEED_BASE};
use crate::rng::Pcg64;

/// Pruning algorithm (`--pruning_algorithm`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algorithm {
    Squared,
    Quartic,
    Wasserstein,
}

impl Algorithm {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "squared" => Some(Algorithm::Squared),
            "quartic" => Some(Algorithm::Quartic),
            "wasserstein" => Some(Algorithm::Wasserstein),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Algorithm::Squared => "squared",
            Algorithm::Quartic => "quartic",
            Algorithm::Wasserstein => "wasserstein",
        }
    }
}

/// Null Z parameters used by the Wasserstein score.
#[derive(Clone, Debug)]
pub struct ZParams {
    pub mean: f64,
    pub scale: f64,
    pub quantiles: Vec<f64>,
    plan_cache: Option<(usize, QuantilePlan)>,
}

impl ZParams {
    pub fn from_null(null: &NullState) -> Self {
        ZParams {
            mean: null.z_mean,
            scale: null.z_scale,
            quantiles: null.z_quantiles.clone(),
            plan_cache: None,
        }
    }

    /// Quantile plan for arrays of `n` values at the interior positions.
    fn plan(&self, n: usize) -> QuantilePlan {
        match &self.plan_cache {
            Some((m, p)) if *m == n => p.clone(),
            _ => QuantilePlan::new(n, &npy::interior_positions(self.quantiles.len())),
        }
    }

    pub fn with_plan_for(mut self, n: usize) -> Self {
        let plan = QuantilePlan::new(n, &npy::interior_positions(self.quantiles.len()));
        self.plan_cache = Some((n, plan));
        self
    }
}

/// The current alignment's score for the algorithm (`initial_global_chi2`).
pub fn initial_score(
    alg: Algorithm,
    expected: &[f64],
    counts: &[i64],
    n_chars: usize,
    nrows: usize,
    z: Option<&ZParams>,
) -> f64 {
    match alg {
        Algorithm::Squared => chi2::global_chi2(expected, counts),
        Algorithm::Quartic => chi2::quartic_global(expected, counts, n_chars, nrows),
        Algorithm::Wasserstein => {
            let z = z.expect("wasserstein needs null Z parameters");
            let rows = chi2::row_chi2(expected, counts, n_chars, nrows);
            let zs: Vec<f64> = rows.iter().map(|&v| (v - z.mean) / z.scale).collect();
            let obs = npy::quantile(&zs, &npy::interior_positions(z.quantiles.len()));
            let diffs: Vec<f64> = obs.iter().zip(&z.quantiles).map(|(o, q)| (o - q).abs()).collect();
            npy::mean(&diffs)
        }
    }
}

/// Block width `B` used by Python for T taxa.
pub fn block_width(nrows: usize) -> usize {
    (2_000_000 / (8 * nrows)).clamp(16, 256).max(16)
}

/// The (start, end) column blocks in Python's chunk/block order.
pub fn blocks(nrows: usize, ncols: usize, num_workers: usize) -> Vec<(usize, usize)> {
    let b = block_width(nrows);
    let n_chunks = if nrows * ncols < 5_000_000 || num_workers <= 1 {
        1
    } else {
        num_workers.min((ncols / b / 4).max(1))
    };
    let chunk_ranges: Vec<(usize, usize)> = if n_chunks > 1 {
        let bounds = npy::linspace_int(ncols, n_chunks + 1);
        bounds.windows(2).map(|w| (w[0], w[1])).collect()
    } else {
        vec![(0, ncols)]
    };
    let mut out = Vec::new();
    for (cs, ce) in chunk_ranges {
        let mut j = cs;
        while j < ce {
            out.push((j, (j + b).min(ce)));
            j += b;
        }
    }
    out
}

/// Shared per-call inputs of the block kernel.
struct Base {
    t_n: usize,
    n_chars: usize,
    cr: Vec<f64>,
    u: Vec<f64>,
    total_chars: Vec<f64>,
    n_total: f64,
    col_totals: Vec<f64>,
}

impl Base {
    fn new(counts: &[i64], t_n: usize, n_chars: usize) -> Self {
        let t = t_n;
        let cr: Vec<f64> = counts.iter().map(|&v| v as f64).collect();
        let u: Vec<f64> = cr.iter().map(|&v| v * v).collect();
        let total_chars: Vec<f64> =
            (0..n_chars).map(|c| npy::sum(&cr[c * t..(c + 1) * t])).collect();
        let n_total = npy::sum(&total_chars);
        let mut col_totals = cr[..t].to_vec();
        for c in 1..n_chars {
            for ti in 0..t {
                col_totals[ti] += cr[c * t + ti];
            }
        }
        Base { t_n, n_chars, cr, u, total_chars, n_total, col_totals }
    }

    /// `_per_row_chi2_block`: per-row chi^2 (T x b, row-major) after removing
    /// each block column (`cols`, the block's b encoded columns in order), plus
    /// the degenerate-column mask.
    fn per_row_block(&self, cols: &[&[i8]]) -> (Vec<f64>, Vec<bool>) {
        let (c_n, t_n, b) = (self.n_chars, self.t_n, cols.len());
        let mut count_at_col = vec![0.0f64; c_n * b];
        let mut n_nongap = vec![0.0f64; b];
        for k in 0..b {
            for &code in cols[k] {
                if code >= 0 {
                    count_at_col[code as usize * b + k] += 1.0;
                    n_nongap[k] += 1.0;
                }
            }
        }
        let n_j: Vec<f64> = n_nongap.iter().map(|&v| self.n_total - v).collect();
        let mut zero_mask = vec![false; b];
        let mut inv_gf = vec![0.0f64; c_n * b];
        for c in 0..c_n {
            for k in 0..b {
                let tc = self.total_chars[c] - count_at_col[c * b + k];
                if tc <= 0.0 {
                    zero_mask[k] = true;
                }
                let safe = if tc > 0.0 { tc } else { 1.0 };
                inv_gf[c * b + k] = n_j[k] / safe;
            }
        }
        let mut a = vec![0.0f64; t_n * b];
        blas::ut_times(&self.u, &inv_gf, c_n, t_n, b, &mut a);
        let mut per_row = vec![0.0f64; t_n * b];
        for k in 0..b {
            let col = cols[k];
            for t in 0..t_n {
                let code = col[t];
                let nongap = if code >= 0 { 1.0 } else { 0.0 };
                let cs = if code >= 0 { code as usize } else { 0 };
                let ct = self.col_totals[t] - nongap;
                if ct <= 0.0 {
                    zero_mask[k] = true;
                }
                let inv_at = inv_gf[cs * b + k];
                let cr_at = self.cr[cs * t_n + t];
                let s = a[t * b + k] + nongap * (1.0 - 2.0 * cr_at) * inv_at;
                let safe_ct = if ct > 0.0 { ct } else { 1.0 };
                per_row[t * b + k] = s / safe_ct - ct;
            }
        }
        (per_row, zero_mask)
    }
}

/// The per-block reduction of the active algorithm.
#[derive(Clone, Copy)]
enum Reducer<'a> {
    Sum,
    SumOfSquares,
    WassersteinZ { z: &'a ZParams, plan: &'a QuantilePlan },
}

/// numpy's `NPY_MIN_ELIDE_BYTES`: temporaries at least this large are reused
/// as the output of the next binary operation (temporary elision).
const MIN_ELIDE_BYTES: usize = 256 * 1024;

/// Whether Python's `per_row_block` for a block of width `b` is F-ordered.
///
/// With an F-ordered `alignment_int`, every intermediate of
/// `_per_row_chi2_block` is F-ordered except the BLAS product `A` (C). In
/// `S = A + tmp`, numpy elides the F-ordered temporary `tmp` (if it is at
/// least 256 KiB) and writes the sum into it, so S and the result stay
/// F-ordered; otherwise the new output takes C order (C wins conflicts).
pub fn per_row_block_is_fortran(alignment_fortran: bool, nrows: usize, b: usize) -> bool {
    alignment_fortran && b >= 2 && nrows * b * 8 >= MIN_ELIDE_BYTES
}

impl<'a> Reducer<'a> {
    /// The reducer of `alg`; Wasserstein needs `z` and a plan for T values.
    fn new(alg: Algorithm, z: Option<&'a ZParams>, plan: Option<&'a QuantilePlan>) -> Self {
        match alg {
            Algorithm::Squared => Reducer::Sum,
            Algorithm::Quartic => Reducer::SumOfSquares,
            Algorithm::Wasserstein => Reducer::WassersteinZ {
                z: z.expect("wasserstein needs null Z parameters"),
                plan: plan.expect("wasserstein needs a quantile plan"),
            },
        }
    }

    /// Reduce a (T x b) row-major block to b values with numpy's operation
    /// order: axis-0 sums are sequential over rows when the block is C-ordered
    /// with b > 1, and pairwise per column when the reduced axis is contiguous
    /// (b == 1, or an F-ordered block). Wasserstein is layout-independent
    /// (np.quantile works on a C-ordered copy).
    fn reduce(&self, per_row: &[f64], t_n: usize, b: usize, fortran: bool) -> Vec<f64> {
        match *self {
            Reducer::Sum | Reducer::SumOfSquares => {
                let sq = matches!(self, Reducer::SumOfSquares);
                let f = |x: f64| if sq { x * x } else { x };
                if b == 1 || fortran {
                    let mut col = vec![0.0f64; t_n];
                    return (0..b)
                        .map(|k| {
                            for t in 0..t_n {
                                col[t] = f(per_row[t * b + k]);
                            }
                            npy::pairwise_sum(&col)
                        })
                        .collect();
                }
                let mut acc: Vec<f64> = per_row[..b].iter().map(|&x| f(x)).collect();
                for t in 1..t_n {
                    let row = &per_row[t * b..(t + 1) * b];
                    for k in 0..b {
                        acc[k] += f(row[k]);
                    }
                }
                acc
            }
            Reducer::WassersteinZ { z, plan } => {
                let kq = z.quantiles.len();
                let mut col = vec![0.0f64; t_n];
                let mut obs = vec![0.0f64; kq];
                let mut diffs = vec![0.0f64; kq];
                // Fast path: x -> (x - mean) / scale is monotone (also in total
                // order, -0.0 before +0.0) when mean is finite and scale is
                // positive and finite, and maps NaN, and only NaN, to NaN. So
                // the sorted z-values at the plan's ranks are the transformed
                // raw values at those ranks: select them on the raw values and
                // transform only those.
                let monotone = z.mean.is_finite() && z.scale.is_finite() && z.scale > 0.0;
                let ranks = plan.ranks();
                let mut at_ranks = vec![0.0f64; ranks.len()];
                let mut scratch = npy::SelectScratch::default();
                (0..b)
                    .map(|k| {
                        if monotone {
                            let mut has_nan = false;
                            for t in 0..t_n {
                                let x = per_row[t * b + k];
                                has_nan |= x.is_nan();
                                col[t] = x;
                            }
                            if has_nan {
                                // np.quantile of an array holding NaN is all NaN.
                                obs.fill(f64::NAN);
                            } else {
                                npy::select_ranks(&col, ranks, &mut at_ranks, &mut scratch);
                                for v in at_ranks.iter_mut() {
                                    *v = (*v - z.mean) / z.scale;
                                }
                                plan.apply_ranked(&at_ranks, &mut obs);
                            }
                        } else {
                            for t in 0..t_n {
                                col[t] = (per_row[t * b + k] - z.mean) / z.scale;
                            }
                            npy::sort_nan_last(&mut col);
                            plan.apply_sorted(&col, &mut obs);
                        }
                        for i in 0..kq {
                            diffs[i] = (obs[i] - z.quantiles[i]).abs();
                        }
                        let total = if b == 1 {
                            npy::pairwise_sum(&diffs)
                        } else {
                            let mut s = diffs[0];
                            for d in &diffs[1..] {
                                s += d;
                            }
                            s
                        };
                        total / kq as f64
                    })
                    .collect()
            }
        }
    }
}

/// Quantile plan for the Wasserstein reducer on T values (`None` otherwise).
fn plan_for(alg: Algorithm, z: Option<&ZParams>, t_n: usize) -> Option<QuantilePlan> {
    match alg {
        Algorithm::Wasserstein => Some(z.expect("wasserstein needs null Z parameters").plan(t_n)),
        _ => None,
    }
}

/// The active score of the alignment without each column of one block
/// (`cols`, in order). Degenerate columns (a character or taxon total that
/// drops to zero) are recomputed directly, as Python does.
fn block_scores(
    base: &Base,
    reducer: &Reducer,
    counts: &[i64],
    cols: &[&[i8]],
    alignment_fortran: bool,
) -> Vec<f64> {
    let (t_n, b) = (base.t_n, cols.len());
    let (per_row, zero) = base.per_row_block(cols);
    let fortran = per_row_block_is_fortran(alignment_fortran, t_n, b);
    let mut vals = reducer.reduce(&per_row, t_n, b, fortran);
    for (k, &is_zero) in zero.iter().enumerate() {
        if is_zero {
            vals[k] = degenerate_score(reducer, counts, base.n_chars, t_n, cols[k]);
        }
    }
    vals
}

/// Direct recomputation of the score without one column.
fn degenerate_score(reducer: &Reducer, counts: &[i64], n_chars: usize, t_n: usize, col: &[i8]) -> f64 {
    let mut cc = counts.to_vec();
    for (ti, &code) in col.iter().enumerate() {
        if code >= 0 {
            cc[code as usize * t_n + ti] -= 1;
        }
    }
    let e = chi2::expected_observed(&cc, n_chars, t_n);
    let per_row = chi2::row_chi2(&e, &cc, n_chars, t_n);
    reducer.reduce(&per_row, t_n, 1, false)[0]
}

/// Delta score per column: `initial - score_without_column_j`
/// (`calculate_{global_chi2,quartic_chi2,wasserstein_zscore}_difference`).
///
/// `counts` are the (fudged) (C, T) counts of the current alignment.
pub fn column_deltas(
    alg: Algorithm,
    counts: &[i64],
    enc: &Encoded,
    n_chars: usize,
    initial: f64,
    z: Option<&ZParams>,
    num_workers_chisq: usize,
    parallel: bool,
) -> Vec<f64> {
    let t_n = enc.nrows;
    let base = Base::new(counts, t_n, n_chars);
    let plan = plan_for(alg, z, t_n);
    let reducer = Reducer::new(alg, z, plan.as_ref());
    let blocks = blocks(t_n, enc.ncols, num_workers_chisq);
    let run = |&(s, e): &(usize, usize)| {
        let cols: Vec<&[i8]> = (s..e).map(|j| enc.column(j)).collect();
        block_scores(&base, &reducer, counts, &cols, enc.fortran)
    };
    let results: Vec<Vec<f64>> = if parallel {
        blocks.par_iter().map(run).collect()
    } else {
        blocks.iter().map(run).collect()
    };
    results.into_iter().flatten().map(|v| initial - v).collect()
}

/// Scores and deltas of the current alignment for the active algorithm
/// (`AlignmentPruner.prune`). Returns (initial score, per-column deltas).
pub fn score_and_deltas(
    alg: Algorithm,
    counts: &[i64],
    expected: &[f64],
    enc: &Encoded,
    n_chars: usize,
    z: Option<&ZParams>,
    num_workers_chisq: usize,
    parallel: bool,
) -> (f64, Vec<f64>) {
    let initial = initial_score(alg, expected, counts, n_chars, enc.nrows, z);
    let deltas = column_deltas(alg, counts, enc, n_chars, initial, z, num_workers_chisq, parallel);
    (initial, deltas)
}

/// Python's `max(deltas.values())` (first maximum; NaN-propagating like `>`).
pub fn python_max(values: &[f64]) -> f64 {
    let mut best = values[0];
    for &v in &values[1..] {
        if v > best {
            best = v;
        }
    }
    best
}

/// Delta-null distribution: for permutation i (seed 12345 + P + i) the max
/// per-column delta of the permuted alignment (`_compute_null_deltas`).
pub fn null_max_deltas(
    alg: Algorithm,
    enc: &Encoded,
    n_chars: usize,
    permutations: usize,
    z: Option<&ZParams>,
) -> Vec<f64> {
    let z = z.map(|zp| zp.clone().with_plan_for(enc.nrows));
    (0..permutations)
        .into_par_iter()
        .map(|i| null_max_delta(alg, enc, n_chars, SEED_BASE + permutations as u64 + i as u64, z.as_ref()))
        .collect()
}

/// Max delta of the alignment permuted with `seed`, without materializing it.
///
/// Pass 1 shuffles column by column and only accumulates the counts. Pass 2
/// restarts the generator from the same seed, so it regenerates exactly the
/// same permuted columns, one kernel block at a time. Python forces a single
/// chunk here (num_workers 1) and the permuted array is C-ordered.
fn null_max_delta(alg: Algorithm, enc: &Encoded, n_chars: usize, seed: u64, z: Option<&ZParams>) -> f64 {
    let t_n = enc.nrows;
    let mut raw = vec![0i64; n_chars * t_n];
    permutation::for_each_permuted_column(enc, seed, |_, col| chi2::add_column_counts(&mut raw, col, t_n));
    let counts = chi2::fudge_counts(raw);
    let initial = {
        let e = chi2::expected_observed(&counts, n_chars, t_n);
        initial_score(alg, &e, &counts, n_chars, t_n, z)
    };
    let base = Base::new(&counts, t_n, n_chars);
    let plan = plan_for(alg, z, t_n);
    let reducer = Reducer::new(alg, z, plan.as_ref());

    let mut rng = Pcg64::from_seed(seed);
    let mut buf = vec![0i8; t_n * block_width(t_n).min(enc.ncols)];
    // Python's max(): first value, then strictly greater (NaN-sticky).
    let mut best: Option<f64> = None;
    for (s, e) in blocks(t_n, enc.ncols, 1) {
        let b = e - s;
        for (k, slot) in buf[..b * t_n].chunks_exact_mut(t_n).enumerate() {
            slot.copy_from_slice(enc.column(s + k));
            rng.shuffle(slot);
        }
        let cols: Vec<&[i8]> = buf[..b * t_n].chunks_exact(t_n).collect();
        for v in block_scores(&base, &reducer, &counts, &cols, false) {
            let d = initial - v;
            best = Some(match best {
                Some(m) if !(d > m) => m,
                _ => d,
            });
        }
    }
    best.expect("alignment has at least one column")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::permutation::{NullMatrix, NullState};

    /// A small encoded alignment with gaps, a column holding the only
    /// occurrences of character 3 (degenerate when removed) and enough columns
    /// for several blocks.
    fn sample(nrows: usize, ncols: usize, n_chars: usize) -> Encoded {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut codes = vec![0i8; nrows * ncols];
        for j in 0..ncols {
            for t in 0..nrows {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let r = (state >> 33) as usize;
                codes[j * nrows + t] = if j == 0 {
                    (n_chars - 1) as i8
                } else if r % 17 == 0 {
                    -1
                } else {
                    ((r + t * (j % 3)) % (n_chars - 1)) as i8
                };
            }
        }
        Encoded::from_codes(codes, nrows, ncols)
    }

    /// The Wasserstein reducer as it was before rank selection: transform
    /// every value, full sort, quantiles of the sorted array.
    fn wasserstein_by_full_sort(per_row: &[f64], t_n: usize, b: usize, z: &ZParams, plan: &QuantilePlan) -> Vec<f64> {
        let kq = z.quantiles.len();
        (0..b)
            .map(|k| {
                let mut col: Vec<f64> = (0..t_n).map(|t| (per_row[t * b + k] - z.mean) / z.scale).collect();
                npy::sort_nan_last(&mut col);
                let mut obs = vec![0.0; kq];
                plan.apply_sorted(&col, &mut obs);
                let diffs: Vec<f64> = obs.iter().zip(&z.quantiles).map(|(o, q)| (o - q).abs()).collect();
                let total = if b == 1 {
                    npy::pairwise_sum(&diffs)
                } else {
                    diffs[1..].iter().fold(diffs[0], |s, d| s + d)
                };
                total / kq as f64
            })
            .collect()
    }

    #[test]
    fn wasserstein_rank_selection_matches_full_sort() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut unit = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for &(t_n, b) in &[(5000usize, 3usize), (12000, 1), (700, 4), (40, 2)] {
            let mut per_row: Vec<f64> = (0..t_n * b)
                .map(|i| {
                    let x = 50.0 + 30.0 * unit();
                    if i % 997 == 0 { x * 1e6 } else { x }
                })
                .collect();
            if b > 1 {
                // A column with NaN, and infinities in another.
                per_row[3 * b + 1] = f64::NAN;
                per_row[5 * b] = f64::INFINITY;
                per_row[9 * b] = f64::NEG_INFINITY;
            }
            let quantiles: Vec<f64> = {
                let mut q: Vec<f64> = (0..200).map(|_| unit() * 4.0 - 2.0).collect();
                q.sort_unstable_by(f64::total_cmp);
                q
            };
            for &(mean, scale) in &[(55.0, 7.5), (65.0, 1e-300), (0.0, 1e300), (55.0, 0.0), (f64::NAN, 1.0)] {
                let z = ZParams { mean, scale, quantiles: quantiles.clone(), plan_cache: None };
                let plan = z.plan(t_n);
                let reducer = Reducer::WassersteinZ { z: &z, plan: &plan };
                let got = reducer.reduce(&per_row, t_n, b, false);
                let want = wasserstein_by_full_sort(&per_row, t_n, b, &z, &plan);
                for (k, (g, w)) in got.iter().zip(&want).enumerate() {
                    assert_eq!(g.to_bits(), w.to_bits(), "T {t_n}, b {b}, col {k}, mean {mean}, scale {scale}");
                }
            }
        }
    }

    #[test]
    fn streamed_delta_null_matches_materialized() {
        let (n_chars, permutations) = (4, 5);
        for &(nrows, ncols) in &[(12usize, 300usize), (40, 97), (1, 20)] {
            let enc = sample(nrows, ncols, n_chars);
            let null = NullState::new(NullMatrix::compute(&enc, n_chars, permutations));
            let z = ZParams::from_null(&null);
            for alg in [Algorithm::Squared, Algorithm::Quartic, Algorithm::Wasserstein] {
                let got = null_max_deltas(alg, &enc, n_chars, permutations, Some(&z));
                let zp = z.clone().with_plan_for(nrows);
                for (i, g) in got.iter().enumerate() {
                    let seed = SEED_BASE + permutations as u64 + i as u64;
                    let p_enc = permutation::permuted(&enc, seed);
                    let counts = chi2::fudge_counts(chi2::count_rows_raw(&p_enc, n_chars));
                    let e = chi2::expected_observed(&counts, n_chars, nrows);
                    let (_, deltas) =
                        score_and_deltas(alg, &counts, &e, &p_enc, n_chars, Some(&zp), 1, false);
                    let want = python_max(&deltas);
                    assert_eq!(g.to_bits(), want.to_bits(), "{alg:?} {nrows}x{ncols} perm {i}");
                }
            }
        }
    }
}
