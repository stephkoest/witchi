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
        QuantilePlan { n, prev, next, gamma }
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
