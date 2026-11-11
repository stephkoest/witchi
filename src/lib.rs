//! Rust port of WitChi v0.2.0 (https://github.com/stephkoest/witchi).
//!
//! The goal is byte-identical output files to the Python implementation; see
//! PLAN.md for the numerical hazards and how each module addresses them.

pub mod alignment;
pub mod blas;
pub mod chi2;
pub mod delta;
pub mod fasta;
pub mod nexus;
pub mod npy;
pub mod permutation;
pub mod pyfmt;
pub mod rng;

/// Version reported in output files (`witchi_version` in the score JSON).
pub const WITCHI_VERSION: &str = "0.2.0";
