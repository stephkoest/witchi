//! Link OpenBLAS for the default `openblas` feature.
//!
//! The library directory comes from `OPENBLAS_LIB_DIR`; by default the
//! `witchi-0.2.0` conda environment's `lib/` is used, so the Rust binary calls
//! the very same OpenBLAS build as Python WitChi (required for bit-identical
//! output). An rpath is embedded so the binary runs without activating conda.

use std::path::PathBuf;

const DEFAULT_LIB_DIR: &str = "/opt/centos7/sw/miniconda3/envs/witchi-0.2.0/lib";

fn main() {
    println!("cargo:rerun-if-env-changed=OPENBLAS_LIB_DIR");
    println!("cargo:rerun-if-changed=build.rs");
    if std::env::var_os("CARGO_FEATURE_OPENBLAS").is_none() {
        return;
    }
    let dir = std::env::var("OPENBLAS_LIB_DIR").unwrap_or_else(|_| DEFAULT_LIB_DIR.to_string());
    let dir = PathBuf::from(dir);
    // conda-forge numpy links libcblas.so.3, a symlink to libopenblasp-*.so;
    // link the same file (it also exports openblas_set_num_threads).
    let lib = std::env::var("OPENBLAS_LIB_NAME").unwrap_or_else(|_| "cblas".to_string());
    println!("cargo:rerun-if-env-changed=OPENBLAS_LIB_NAME");
    if !dir.join(format!("lib{lib}.so")).exists() {
        panic!(
            "lib{lib}.so not found in {} — set OPENBLAS_LIB_DIR / OPENBLAS_LIB_NAME, or build \
             with `--no-default-features` for the pure-Rust BLAS fallback",
            dir.display()
        );
    }
    println!("cargo:rustc-link-search=native={}", dir.display());
    println!("cargo:rustc-link-lib=dylib={lib}");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
}
