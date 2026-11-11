//! Rust port of WitChi v0.2.0 (https://github.com/stephkoest/witchi).
//!
//! The command-line interface mirrors the Python `argparse` definition in
//! `witchi/witchi.py`, including its underscore-style flag names.

use anyhow::{Result, bail};
use clap::{Args, Parser, Subcommand};

/// Version reported in output files (`witchi_version` in the score JSON).
pub const WITCHI_VERSION: &str = "0.2.0";

#[derive(Parser, Debug)]
#[command(
    name = "witchi",
    version = WITCHI_VERSION,
    about = "Recursively prune alignment based on chi-square test."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Prune alignment based on chi-square test.
    Prune(PruneArgs),
    /// Run permutation test on alignment.
    Test(TestArgs),
}

#[derive(Args, Debug)]
struct PruneArgs {
    /// File containing the alignment
    #[arg(long)]
    file: String,

    /// Alignment format, default is fasta
    #[arg(long, default_value = "fasta")]
    format: String,

    /// Maximum number of columns to prune
    #[arg(long = "max_residue_pruned", default_value_t = 100)]
    max_residue_pruned: i64,

    /// Number of alignment permutations.
    #[arg(long, default_value_t = 100)]
    permutations: usize,

    /// Number of parallel workers (cores) to use for column delta score calculation.
    #[arg(long = "num_workers_chisq", default_value_t = 1)]
    num_workers_chisq: usize,

    /// Number of parallel workers (cores) to use for permutation testing.
    #[arg(long = "num_workers_permute", default_value_t = 1)]
    num_workers_permute: usize,

    /// Number of top delta-chi2 columns to remove.
    #[arg(long = "top_n", default_value_t = 1)]
    top_n: usize,

    /// Pruning algorithm to use: squared, quartic, wasserstein
    #[arg(long = "pruning_algorithm", default_value = "wasserstein")]
    pruning_algorithm: String,

    /// Test per-column delta significance against a permutation null (default).
    #[arg(long = "delta-null", overrides_with = "no_delta_null")]
    delta_null: bool,

    /// Disable the delta-null stopping criterion.
    #[arg(long = "no-delta-null", overrides_with = "delta_null")]
    no_delta_null: bool,
}

#[derive(Args, Debug)]
struct TestArgs {
    /// File containing the alignment
    #[arg(long)]
    file: String,

    /// Alignment format. (fasta)
    #[arg(long, default_value = "fasta")]
    format: String,

    /// Number of parallel workers (cores) to use. (1)
    #[arg(long = "num_workers_permute", default_value_t = 1)]
    num_workers_permute: usize,

    /// Number of alignment permutations. (100)
    #[arg(long, default_value_t = 100)]
    permutations: usize,

    /// Flag to create output file with scores.
    #[arg(long = "create_output")]
    create_output: bool,
}

impl PruneArgs {
    /// argparse.BooleanOptionalAction semantics: last flag wins, default on.
    fn delta_null_enabled(&self) -> bool {
        !self.no_delta_null
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    println!("witchi {WITCHI_VERSION}");
    match cli.command {
        Command::Prune(args) => {
            let _ = args.delta_null_enabled();
            bail!("prune is not implemented yet (args: {args:?})")
        }
        Command::Test(args) => bail!("test is not implemented yet (args: {args:?})"),
    }
}
