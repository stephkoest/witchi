//! Permutation null distribution (`PermutationTest` in Python) and the
//! derived null state used by the pruner.

use rayon::prelude::*;

use crate::chi2::{self, Encoded};
use crate::npy;
use crate::rng::Pcg64;

/// Base seed used by WitChi: permutation `i` uses `12345 + i`.
pub const SEED_BASE: u64 = 12345;

/// Visit the columns of a column-major encoded alignment after shuffling each
/// one with a single generator seeded by `seed`, in column order
/// (`np.apply_along_axis(default_rng(seed).permutation, 0, a)`).
pub fn for_each_permuted_column(enc: &Encoded, seed: u64, mut visit: impl FnMut(usize, &[i8])) {
    let mut rng = Pcg64::from_seed(seed);
    let mut col = vec![0i8; enc.nrows];
    for j in 0..enc.ncols {
        col.copy_from_slice(enc.column(j));
        rng.shuffle(&mut col);
        visit(j, &col);
    }
}

/// The whole permuted alignment (column-major), for the delta null.
pub fn permuted(enc: &Encoded, seed: u64) -> Encoded {
    let mut codes = vec![0i8; enc.codes.len()];
    let t = enc.nrows;
    for_each_permuted_column(enc, seed, |j, col| codes[j * t..(j + 1) * t].copy_from_slice(col));
    Encoded { nrows: enc.nrows, ncols: enc.ncols, codes, fortran: false }
}

/// Per-taxon chi^2 of one permuted alignment.
pub fn permuted_row_chi2(enc: &Encoded, n_chars: usize, seed: u64) -> Vec<f64> {
    let mut counts = vec![0i64; n_chars * enc.nrows];
    for_each_permuted_column(enc, seed, |_, col| chi2::add_column_counts(&mut counts, col, enc.nrows));
    let counts = chi2::fudge_counts(counts);
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
    pub fn compute(enc: &Encoded, n_chars: usize, permutations: usize, threads: usize) -> Self {
        let rows: Vec<Vec<f64>> = with_threads(threads, || {
            (0..permutations)
                .into_par_iter()
                .map(|i| permuted_row_chi2(enc, n_chars, SEED_BASE + i as u64))
                .collect()
        });
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

/// Per-taxon Bonferroni empirical p-values, clipped to [0, 1].
pub fn empirical_p_per_taxon(per_row: &[f64], null: &[f64]) -> Vec<f64> {
    let n_taxa = per_row.len() as f64;
    per_row
        .iter()
        .map(|&x| {
            let p = (null.iter().filter(|&&v| v >= x).count() as f64 / null.len() as f64) * n_taxa;
            p.max(0.0).min(1.0)
        })
        .collect()
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
