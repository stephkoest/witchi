//! Chi-squared building blocks (`ChiSquareCalculator` in Python).
//!
//! Layout conventions (matching numpy so reduction order can be read off):
//! * alignments are encoded column-major: `codes[j * nrows + t]`, value =
//!   index into the char set or -1 (gap / not in the char set);
//! * counts and expected values are (C, T) row-major: `x[c * T + t]`.

use crate::npy;

/// An int-encoded alignment, column-major.
#[derive(Clone, Debug)]
pub struct Encoded {
    pub nrows: usize,
    pub ncols: usize,
    pub codes: Vec<i8>,
}

impl Encoded {
    /// Encode a row-major byte alignment against `char_set`.
    pub fn from_rows(data: &[u8], nrows: usize, ncols: usize, char_set: &[u8]) -> Self {
        let lut = lookup_table(char_set);
        let mut codes = vec![-1i8; nrows * ncols];
        for t in 0..nrows {
            let row = &data[t * ncols..(t + 1) * ncols];
            for (j, &b) in row.iter().enumerate() {
                codes[j * nrows + t] = lut[b as usize];
            }
        }
        Encoded { nrows, ncols, codes }
    }

    pub fn column(&self, j: usize) -> &[i8] {
        &self.codes[j * self.nrows..(j + 1) * self.nrows]
    }

    /// Remove the given columns (any order; duplicates ignored).
    pub fn delete_columns(&self, cols: &[usize]) -> Encoded {
        let mut drop = vec![false; self.ncols];
        for &c in cols {
            drop[c] = true;
        }
        let keep = drop.iter().filter(|d| !**d).count();
        let mut codes = Vec::with_capacity(keep * self.nrows);
        for j in 0..self.ncols {
            if !drop[j] {
                codes.extend_from_slice(self.column(j));
            }
        }
        Encoded { nrows: self.nrows, ncols: keep, codes }
    }
}

/// Byte -> char-set index (or -1).
pub fn lookup_table(char_set: &[u8]) -> [i8; 256] {
    let mut lut = [-1i8; 256];
    for (i, &c) in char_set.iter().enumerate() {
        lut[c as usize] = i as i8;
    }
    lut
}

/// Add the residue counts of one encoded column into (C, T) `counts`.
#[inline]
pub fn add_column_counts(counts: &mut [i64], col: &[i8], nrows: usize) {
    for (t, &c) in col.iter().enumerate() {
        if c >= 0 {
            counts[c as usize * nrows + t] += 1;
        }
    }
}

/// Un-fudged (C, T) counts (`_count_rows_from_int`).
pub fn count_rows_raw(enc: &Encoded, n_chars: usize) -> Vec<i64> {
    let mut counts = vec![0i64; n_chars * enc.nrows];
    for j in 0..enc.ncols {
        add_column_counts(&mut counts, enc.column(j), enc.nrows);
    }
    counts
}

/// The +1 zero-safety fudge of `calculate_row_counts` / `_cached_count_rows`.
pub fn fudge_counts(mut counts: Vec<i64>) -> Vec<i64> {
    if counts.iter().any(|&c| c == 0) {
        for c in counts.iter_mut() {
            *c += 1;
        }
    }
    counts
}

/// `calculate_expected_observed`: e[c,t] = gf[c] * col_total[t] (+1 fudge).
pub fn expected_observed(counts: &[i64], n_chars: usize, nrows: usize) -> Vec<f64> {
    let total_chars: Vec<i64> =
        (0..n_chars).map(|c| counts[c * nrows..(c + 1) * nrows].iter().sum()).collect();
    let total: i64 = total_chars.iter().sum();
    let gf: Vec<f64> = total_chars.iter().map(|&tc| tc as f64 / total as f64).collect();
    let mut col_totals = vec![0i64; nrows];
    for c in 0..n_chars {
        for t in 0..nrows {
            col_totals[t] += counts[c * nrows + t];
        }
    }
    let mut e = vec![0.0f64; n_chars * nrows];
    let mut any_zero = false;
    for c in 0..n_chars {
        for t in 0..nrows {
            let v = gf[c] * col_totals[t] as f64;
            any_zero |= v == 0.0;
            e[c * nrows + t] = v;
        }
    }
    if any_zero {
        for v in e.iter_mut() {
            *v += 1.0;
        }
    }
    e
}

/// Elementwise (o - e)^2 / e in (C, T) layout.
pub fn chi2_values(expected: &[f64], counts: &[i64]) -> Vec<f64> {
    counts
        .iter()
        .zip(expected)
        .map(|(&o, &e)| {
            let d = o as f64 - e;
            d * d / e
        })
        .collect()
}

/// `calculate_row_chi2`: per-taxon chi^2, summed over characters sequentially
/// (axis-0 reduction of a C-contiguous (C, T) array).
pub fn row_chi2(expected: &[f64], counts: &[i64], n_chars: usize, nrows: usize) -> Vec<f64> {
    let vals = chi2_values(expected, counts);
    let mut acc = vals[..nrows].to_vec();
    for c in 1..n_chars {
        for t in 0..nrows {
            acc[t] += vals[c * nrows + t];
        }
    }
    if nrows == 1 {
        // (C, 1): numpy reduces the contiguous length-C axis pairwise.
        acc[0] = npy::pairwise_sum(&vals);
    }
    acc
}

/// `calculate_global_chi2`: sum of all (C, T) chi^2 values (pairwise).
pub fn global_chi2(expected: &[f64], counts: &[i64]) -> f64 {
    npy::sum(&chi2_values(expected, counts))
}

/// `calculate_quartic_row_global_chi2`: sum over taxa of chi^2 squared.
pub fn quartic_global(expected: &[f64], counts: &[i64], n_chars: usize, nrows: usize) -> f64 {
    let rows = row_chi2(expected, counts, n_chars, nrows);
    let sq: Vec<f64> = rows.iter().map(|x| x * x).collect();
    npy::sum(&sq)
}

/// Counts with the fudge, expected values and per-taxon chi^2 of an encoded
/// alignment (`calculate_row_counts` → `calculate_expected_observed` →
/// `calculate_row_chi2`).
pub fn per_row_chi2_of(enc: &Encoded, n_chars: usize) -> Vec<f64> {
    let counts = fudge_counts(count_rows_raw(enc, n_chars));
    let e = expected_observed(&counts, n_chars, enc.nrows);
    row_chi2(&e, &counts, n_chars, enc.nrows)
}
