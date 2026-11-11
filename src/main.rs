//! Rust port of WitChi v0.2.0 (https://github.com/stephkoest/witchi).
//!
//! The command-line interface mirrors the Python `argparse` definition in
//! `witchi/witchi.py`, including its underscore-style flag names.

use anyhow::Result;
use clap::{Args, Parser, Subcommand};


use witchi::WITCHI_VERSION;

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

    /// Rust port only: number of worker threads. Defaults to the larger of
    /// --num_workers_chisq and --num_workers_permute. Does not affect results
    /// (unlike --num_workers_chisq, which sets Python's column chunking).
    #[arg(long, default_value_t = 0)]
    threads: usize,
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

    /// Rust port only: number of worker threads (default: --num_workers_permute).
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

impl PruneArgs {
    /// argparse.BooleanOptionalAction semantics: last flag wins, default on.
    fn delta_null_enabled(&self) -> bool {
        !self.no_delta_null
    }
}

/// OpenBLAS starts its own worker threads when the library is loaded, which
/// then spin and waste CPU even though every call here is single-threaded.
/// It only honours `OPENBLAS_NUM_THREADS` at load time, so re-exec once with
/// it set (results are unaffected: the same single-threaded kernels run).
#[cfg(all(feature = "openblas", unix))]
fn ensure_single_threaded_openblas() {
    use std::os::unix::process::CommandExt;
    if std::env::var_os("OPENBLAS_NUM_THREADS").is_some() {
        return;
    }
    if let Ok(exe) = std::env::current_exe() {
        let err = std::process::Command::new(exe)
            .args(std::env::args_os().skip(1))
            .env("OPENBLAS_NUM_THREADS", "1")
            .exec();
        eprintln!("warning: could not re-exec with OPENBLAS_NUM_THREADS=1: {err}");
    }
}

fn main() -> Result<()> {
    #[cfg(all(feature = "openblas", unix))]
    ensure_single_threaded_openblas();
    let cli = Cli::parse();
    println!("witchi {WITCHI_VERSION}");
    match cli.command {
        Command::Prune(args) => {
            let opts = witchi::pruner::PruneOptions {
                delta_null: args.delta_null_enabled(),
                file: args.file,
                format: args.format,
                max_residue_pruned: args.max_residue_pruned,
                permutations: args.permutations,
                num_workers_chisq: args.num_workers_chisq,
                num_workers_permute: args.num_workers_permute,
                top_n: args.top_n,
                pruning_algorithm: args.pruning_algorithm,
                threads: args.threads,
            };
            witchi::pruner::run(&opts)
        }
        Command::Test(args) => witchi::testcmd::run(&witchi::testcmd::TestOptions {
            file: args.file,
            format: args.format,
            num_workers_permute: args.num_workers_permute,
            permutations: args.permutations,
            create_output: args.create_output,
            threads: args.threads,
        }),
    }
}
