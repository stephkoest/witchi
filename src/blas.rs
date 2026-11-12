//! The one BLAS product WitChi computes: `U.T @ inv_gf`, with U (C, T) and
//! inv_gf (C, b) both C-contiguous, giving (T, b) C-contiguous.
//!
//! numpy 2.4.6 (`umath/matmul.c.src`) dispatches this as:
//! * T > 1, b > 1: `cblas_dgemm(RowMajor, Trans, NoTrans, T, b, C, 1, U, T, inv_gf, b, 0, out, b)`
//! * T > 1, b == 1: `cblas_dgemv(RowMajor, Trans, C, T, 1, U, T, inv_gf, 1, 0, out, 1)`
//! * T == 1, b > 1: `cblas_dgemv(RowMajor, Trans, C, b, 1, inv_gf, b, U, 1, 0, out, 1)`
//! * T == 1, b == 1: `0.0 + cblas_ddot(C, U, 1, inv_gf, 1)`
//!
//! The `openblas` feature (default) calls OpenBLAS with exactly these
//! arguments, so results are bit-identical to Python on the same machine and
//! OpenBLAS build. Without it, a portable FMA chain is used (matches OpenBLAS's
//! interior gemm tiles on AVX2/FMA machines, but not every edge case).

/// out (T x b, row-major) = U^T (T x C) * inv_gf (C x b).
pub fn ut_times(u: &[f64], inv_gf: &[f64], c: usize, t: usize, b: usize, out: &mut [f64]) {
    debug_assert_eq!(u.len(), c * t);
    debug_assert_eq!(inv_gf.len(), c * b);
    debug_assert_eq!(out.len(), t * b);
    backend::ut_times(u, inv_gf, c, t, b, out)
}

/// Make the BLAS library single-threaded (Python wraps the kernel in
/// `threadpool_limits(limits=1, user_api="blas")`).
pub fn init_single_threaded() {
    backend::init();
}

/// Name of the compiled-in backend (for diagnostics).
pub fn backend_name() -> &'static str {
    backend::NAME
}

#[cfg(feature = "openblas")]
mod backend {
    use std::os::raw::c_int;

    const ROW_MAJOR: c_int = 101;
    const NO_TRANS: c_int = 111;
    const TRANS: c_int = 112;

    unsafe extern "C" {
        fn cblas_dgemm(
            order: c_int, transa: c_int, transb: c_int, m: c_int, n: c_int, k: c_int,
            alpha: f64, a: *const f64, lda: c_int, b: *const f64, ldb: c_int, beta: f64,
            c: *mut f64, ldc: c_int,
        );
        fn cblas_dgemv(
            order: c_int, trans: c_int, m: c_int, n: c_int, alpha: f64, a: *const f64,
            lda: c_int, x: *const f64, incx: c_int, beta: f64, y: *mut f64, incy: c_int,
        );
        fn cblas_ddot(n: c_int, x: *const f64, incx: c_int, y: *const f64, incy: c_int) -> f64;
        fn openblas_set_num_threads(n: c_int);
    }

    pub const NAME: &str = "openblas";

    pub fn init() {
        unsafe { openblas_set_num_threads(1) };
    }

    pub fn ut_times(u: &[f64], inv_gf: &[f64], c: usize, t: usize, b: usize, out: &mut [f64]) {
        let (ci, ti, bi) = (c as c_int, t as c_int, b as c_int);
        unsafe {
            if t > 1 && b > 1 {
                cblas_dgemm(
                    ROW_MAJOR, TRANS, NO_TRANS, ti, bi, ci, 1.0, u.as_ptr(), ti,
                    inv_gf.as_ptr(), bi, 0.0, out.as_mut_ptr(), bi,
                );
            } else if b == 1 && t > 1 {
                cblas_dgemv(
                    ROW_MAJOR, TRANS, ci, ti, 1.0, u.as_ptr(), ti, inv_gf.as_ptr(), 1, 0.0,
                    out.as_mut_ptr(), 1,
                );
            } else if t == 1 && b > 1 {
                cblas_dgemv(
                    ROW_MAJOR, TRANS, ci, bi, 1.0, inv_gf.as_ptr(), bi, u.as_ptr(), 1, 0.0,
                    out.as_mut_ptr(), 1,
                );
            } else {
                out[0] = 0.0 + cblas_ddot(ci, u.as_ptr(), 1, inv_gf.as_ptr(), 1);
            }
        }
    }
}

#[cfg(not(feature = "openblas"))]
mod backend {
    pub const NAME: &str = "pure-rust";

    pub fn init() {}

    pub fn ut_times(u: &[f64], inv_gf: &[f64], c: usize, t: usize, b: usize, out: &mut [f64]) {
        for ti in 0..t {
            for k in 0..b {
                let mut acc = 0.0f64;
                for ci in 0..c {
                    acc = u[ci * t + ti].mul_add(inv_gf[ci * b + k], acc);
                }
                out[ti * b + k] = acc;
            }
        }
    }
}
