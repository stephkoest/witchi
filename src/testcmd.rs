//! `witchi test`: permutation test of an alignment
//! (`PermutationTest.run_test` in Python).

use anyhow::{Result, bail};
use indexmap::IndexMap;
use std::path::Path;
use std::time::Instant;

use crate::alignment;
use crate::chi2::{self, Encoded};
use crate::npy;
use crate::permutation::{self, NullMatrix, RobustStats};
use crate::pyfmt::{fixed, repr_f64, tsv_row};

/// Command-line options of `witchi test`.
#[derive(Clone, Debug)]
pub struct TestOptions {
    pub file: String,
    pub format: String,
    pub num_workers_permute: usize,
    pub permutations: usize,
    pub create_output: bool,
    /// Rust-only: worker threads (0 = `num_workers_permute`).
    pub threads: usize,
}

/// Python's output name: `file.replace("." + file.split(".")[-1], "_scores.tsv")`
/// (note: replaces *every* occurrence of the final extension).
pub fn scores_path(file: &str) -> String {
    let ext = file.rsplit('.').next().unwrap_or(file);
    file.replace(&format!(".{ext}"), "_scores.tsv")
}

pub fn run(opts: &TestOptions) -> Result<()> {
    let threads = if opts.threads > 0 { opts.threads } else { opts.num_workers_permute.max(1) };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| anyhow::anyhow!("cannot start thread pool: {e}"))?;
    pool.install(|| run_inner(opts))
}

fn run_inner(opts: &TestOptions) -> Result<()> {
    let start = Instant::now();
    if opts.permutations == 0 {
        bail!("--permutations must be at least 1");
    }
    let out_path = scores_path(&opts.file);
    if opts.create_output && out_path == opts.file {
        // Python would overwrite the input alignment with the scores table.
        bail!(
            "cannot derive an output name from {:?} (no extension); Python WitChi would overwrite the input",
            opts.file
        );
    }
    let aln = alignment::read_alignment(Path::new(&opts.file), &opts.format)?;
    let seq_type = alignment::detect_sequence_type(&aln)?;
    println!("{}", seq_type.message());
    let char_set = seq_type.char_set();
    let n_chars = char_set.len();
    let enc = Encoded::from_rows(&aln.data, aln.nrows(), aln.ncols, char_set);

    println!("Running {} permutations.", opts.permutations);
    println!("Using {} worker(s) for permutation", opts.num_workers_permute);
    let null = NullMatrix::compute(&enc, n_chars, opts.permutations);
    let sums = &null.sums;
    let pooled = &null.pooled;

    let per_row = chi2::per_row_chi2_of(&enc, n_chars);
    let total = npy::sum(&per_row);
    let alignment_z = RobustStats::of(sums).zscore(total);
    let min_sum = sums.iter().copied().fold(f64::NAN, |m, v| if m.is_nan() || v < m { v } else { m });
    let max_sum = sums.iter().copied().fold(f64::NAN, |m, v| if m.is_nan() || v > m { v } else { m });
    println!(
        "Alignment chi2score: {} | Permutations alignment chi2scores: {} - {} | Empirical-P: {}",
        fixed(total, 2),
        fixed(min_sum, 2),
        fixed(max_sum, 2),
        repr_f64(permutation::empirical_p(total, sums))
    );
    println!(
        "Alignment mean taxa chi2score: {} | Permutations mean taxa chi2score: {} | Alignment Z-score: {}",
        fixed(npy::mean(&per_row), 2),
        fixed(npy::mean(pooled), 2),
        fixed(alignment_z, 2)
    );

    // make_score_dict: {taxon: (p, z)} (dict semantics for duplicate ids),
    // sorted by z descending (stable).
    let pvals = permutation::empirical_p_per_taxon(&per_row, pooled);
    let rs = RobustStats::of(pooled);
    let mut by_name: IndexMap<&str, (f64, f64)> = IndexMap::new();
    for (i, rec) in aln.records.iter().enumerate() {
        by_name.insert(rec.id.as_str(), (pvals[i], rs.zscore(per_row[i])));
    }
    let mut scores: Vec<(&str, f64, f64)> = by_name.into_iter().map(|(k, (p, z))| (k, p, z)).collect();
    scores.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    let significant = scores.iter().filter(|s| s.1 < 0.05).count();
    println!(
        "Biased Taxa based on corrected empirical-P-values: {significant} of {}",
        per_row.len()
    );
    if opts.create_output {
        let mut by_abs = scores.clone();
        by_abs.sort_by(|a, b| b.2.abs().partial_cmp(&a.2.abs()).unwrap_or(std::cmp::Ordering::Equal));
        let mut tsv = String::new();
        tsv_row(&mut tsv, &["Row", "Empirical-Pvalue", "Z-Score"]);
        for (name, p, z) in &by_abs {
            tsv_row(&mut tsv, &[name, &fixed(*p, 4), &fixed(*z, 2)]);
        }
        std::fs::write(&out_path, tsv.as_bytes())?;
        println!("Printing taxa p-values and z-scores to file: {out_path}");
        println!("Top 5 taxa z-scores and corresponding empirical p-values:");
        for (name, p, z) in scores.iter().take(5) {
            println!("{name}\t{}\t{}", fixed(*z, 2), fixed(*p, 4));
        }
    } else {
        println!("Taxa z-scores and empirical p-values:");
        for (name, p, z) in &scores {
            println!("{name}\t{}\t{}", fixed(*z, 2), fixed(*p, 4));
        }
    }
    println!("Execution time for testing: {} seconds", fixed(start.elapsed().as_secs_f64(), 2));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::scores_path;

    #[test]
    fn python_replace_semantics() {
        assert_eq!(scores_path("aln.fasta"), "aln_scores.tsv");
        assert_eq!(scores_path("dir/aln.faa"), "dir/aln_scores.tsv");
        assert_eq!(scores_path("x.faa/aln.faa"), "x_scores.tsv/aln_scores.tsv");
        assert_eq!(scores_path("noext"), "noext");
    }
}
