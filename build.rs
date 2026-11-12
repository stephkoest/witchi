//! Link OpenBLAS for the default `openblas` feature.
//!
//! The library directory is `OPENBLAS_LIB_DIR` if set, otherwise the active
//! conda environment's `lib/` (`$CONDA_PREFIX/lib`), e.g. the environment
//! created from `environment.yml`. Linking the same OpenBLAS build that numpy
//! uses is what makes the output bit-identical to Python WitChi. An rpath is
//! embedded so the binary runs without activating the environment.

use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-env-changed=OPENBLAS_LIB_DIR");
    println!("cargo:rerun-if-env-changed=OPENBLAS_LIB_NAME");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    println!("cargo:rerun-if-changed=build.rs");
    if std::env::var_os("CARGO_FEATURE_OPENBLAS").is_none() {
        return;
    }
    let dir = match std::env::var_os("OPENBLAS_LIB_DIR") {
        Some(d) => PathBuf::from(d),
        None => match std::env::var_os("CONDA_PREFIX") {
            Some(p) => PathBuf::from(p).join("lib"),
            None => panic!(
                "no OpenBLAS location: activate the conda environment from environment.yml, \
                 set OPENBLAS_LIB_DIR, or build with `--no-default-features` for the \
                 pure-Rust BLAS fallback"
            ),
        },
    };
    // conda-forge numpy links libcblas.so.3, a symlink to libopenblasp-*.so;
    // link the same file (it also exports openblas_set_num_threads).
    let lib = std::env::var("OPENBLAS_LIB_NAME").unwrap_or_else(|_| "cblas".to_string());
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
