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
struct Base<'a> {
    enc: &'a Encoded,
    n_chars: usize,
    cr: Vec<f64>,
    u: Vec<f64>,
    total_chars: Vec<f64>,
    n_total: f64,
    col_totals: Vec<f64>,
}

impl<'a> Base<'a> {
    fn new(counts: &[i64], enc: &'a Encoded, n_chars: usize) -> Self {
        let t = enc.nrows;
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
        Base { enc, n_chars, cr, u, total_chars, n_total, col_totals }
    }

    /// `_per_row_chi2_block`: per-row chi^2 (T x b, row-major) after removing
    /// each block column, plus the degenerate-column mask.
    fn per_row_block(&self, j_start: usize, j_end: usize) -> (Vec<f64>, Vec<bool>) {
        let (c_n, t_n, b) = (self.n_chars, self.enc.nrows, j_end - j_start);
        let mut count_at_col = vec![0.0f64; c_n * b];
        let mut n_nongap = vec![0.0f64; b];
        for k in 0..b {
            for &code in self.enc.column(j_start + k) {
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
            let col = self.enc.column(j_start + k);
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

impl Reducer<'_> {
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
                (0..b)
                    .map(|k| {
                        for t in 0..t_n {
                            col[t] = (per_row[t * b + k] - z.mean) / z.scale;
                        }
                        npy::sort_nan_last(&mut col);
                        plan.apply_sorted(&col, &mut obs);
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
    let base = Base::new(counts, enc, n_chars);
    let t_n = enc.nrows;
    let plan;
    let reducer = match alg {
        Algorithm::Squared => Reducer::Sum,
        Algorithm::Quartic => Reducer::SumOfSquares,
        Algorithm::Wasserstein => {
            let zp = z.expect("wasserstein needs null Z parameters");
            plan = zp.plan(t_n);
            Reducer::WassersteinZ { z: zp, plan: &plan }
        }
    };
    let blocks = blocks(t_n, enc.ncols, num_workers_chisq);
    let run = |&(s, e): &(usize, usize)| {
        let (per_row, zero) = base.per_row_block(s, e);
        let fortran = per_row_block_is_fortran(enc.fortran, t_n, e - s);
        (reducer.reduce(&per_row, t_n, e - s, fortran), zero)
    };
    let results: Vec<(Vec<f64>, Vec<bool>)> = if parallel {
        blocks.par_iter().map(run).collect()
    } else {
        blocks.iter().map(run).collect()
    };
    let mut out = Vec::with_capacity(enc.ncols);
    let mut zero_cols = Vec::new();
    for ((s, _), (vals, zero)) in blocks.iter().zip(results) {
        for (k, &is_zero) in zero.iter().enumerate() {
            if is_zero {
                zero_cols.push(s + k);
            }
        }
        out.extend(vals);
    }
    // Degenerate columns: direct recomputation without column j.
    for j in zero_cols {
        let mut cc = counts.to_vec();
        for (ti, &code) in enc.column(j).iter().enumerate() {
            if code >= 0 {
                cc[code as usize * t_n + ti] -= 1;
            }
        }
        let e = chi2::expected_observed(&cc, n_chars, t_n);
        let per_row = chi2::row_chi2(&e, &cc, n_chars, t_n);
        out[j] = reducer.reduce(&per_row, t_n, 1, false)[0];
    }
    out.iter().map(|&v| initial - v).collect()
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
    threads: usize,
) -> Vec<f64> {
    let z = z.map(|zp| zp.clone().with_plan_for(enc.nrows));
    permutation::with_threads(threads, || {
        (0..permutations)
            .into_par_iter()
            .map(|i| {
                let seed = SEED_BASE + permutations as u64 + i as u64;
                let p_enc = permutation::permuted(enc, seed);
                let counts = chi2::fudge_counts(chi2::count_rows_raw(&p_enc, n_chars));
                let e = chi2::expected_observed(&counts, n_chars, enc.nrows);
                // num_workers is forced to 1 in Python: a single chunk.
                let (_, deltas) =
                    score_and_deltas(alg, &counts, &e, &p_enc, n_chars, z.as_ref(), 1, false);
                python_max(&deltas)
            })
            .collect()
    })
}
