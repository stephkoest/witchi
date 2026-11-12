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
pub mod output;
pub mod permutation;
pub mod pruner;
pub mod pyfmt;
pub mod rng;
pub mod testcmd;

/// Version reported in the console banner and the JSON `witchi_version` field
/// (the crate version; the Python reference these outputs match is 0.2.0).
pub const WITCHI_VERSION: &str = env!("CARGO_PKG_VERSION");
