//! The pruning workflow (`AlignmentPruner` in Python): a literal port of
//! `run` and `recursive_prune`, including touchdown rollback, the one-shot
//! null recompute, delta-null walk and the cap bookkeeping.

use anyhow::{Result, bail};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use crate::alignment::{self, Alignment, SequenceType};
use crate::chi2::{self, Encoded};
use crate::delta::{self, Algorithm, ZParams};
use crate::output;
use crate::permutation::{self, NullMatrix, NullState, RobustStats};
use crate::pyfmt::fixed;

/// Command-line options of `witchi prune`.
#[derive(Clone, Debug)]
pub struct PruneOptions {
    pub file: String,
    pub format: String,
    pub max_residue_pruned: i64,
    pub permutations: usize,
    pub num_workers_chisq: usize,
    pub num_workers_permute: usize,
    pub top_n: usize,
    pub pruning_algorithm: String,
    pub delta_null: bool,
}

/// One row of the prune trace (`prune_dict` value).
#[derive(Clone, Debug)]
pub struct PruneRow {
    pub iteration: usize,
    pub original_index: usize,
    pub global_chi2: f64,
    pub initial_score: f64,
    pub difference: f64,
    pub significant_taxa: usize,
    pub null_pvalue: Option<f64>,
}

/// Alignment-level summary against the original null (`_reported_alignment_stats`).
#[derive(Clone, Copy, Debug)]
pub struct ReportedStats {
    pub z: f64,
    pub p: f64,
    pub biased: usize,
}

/// Everything `run` writes out.
pub struct PruneResult {
    /// Original column indices kept, in order.
    pub kept_columns: Vec<usize>,
    pub prune_dict: IndexMap<usize, PruneRow>,
    pub before_real: Vec<f64>,
    pub after_real: Vec<f64>,
    pub stop_reason: String,
    pub null_recomputes: usize,
    pub null_max_deltas: Option<Vec<f64>>,
    pub before_stats: ReportedStats,
    pub after_stats: ReportedStats,
}

/// Current null plus the delta-null derived from it.
struct OperativeNull {
    state: NullState,
    z: ZParams,
    max_deltas: Option<Vec<f64>>,
}

struct Pruner<'a> {
    opts: &'a PruneOptions,
    alg: Algorithm,
    n_chars: usize,
    top_n: usize,
    initial_top_n: usize,
    null_recompute_count: usize,
    last_recompute_ncols: Option<usize>,
}

/// Per-iteration statistics of the current alignment (`_calculate_per_row_stats`).
struct Stats {
    counts: Vec<i64>,
    expected: Vec<f64>,
    per_row: Vec<f64>,
    sum: f64,
}

impl Pruner<'_> {
    fn compute_null_state(&self, enc: &Encoded) -> OperativeNull {
        let null = NullMatrix::compute(enc, self.n_chars, self.opts.permutations, self.opts.num_workers_permute);
        let state = NullState::new(null);
        let z = ZParams::from_null(&state);
        let max_deltas = if self.opts.delta_null {
            println!("Computing null delta distribution for stopping criterion...");
            let md = delta::null_max_deltas(
                self.alg,
                enc,
                self.n_chars,
                self.opts.permutations,
                Some(&z),
                self.opts.num_workers_permute,
            );
            println!(
                "Null max-delta distribution (n={}): mean={}, p95={}",
                self.opts.permutations,
                fixed(crate::npy::mean(&md), 6),
                fixed(crate::npy::percentile(&md, 95.0), 6)
            );
            Some(md)
        } else {
            None
        };
        OperativeNull { state, z, max_deltas }
    }

    fn recompute_due(&self, ncols: usize, iteration: usize) -> bool {
        if iteration == 0 || self.null_recompute_count >= 1 {
            return false;
        }
        Some(ncols) != self.last_recompute_ncols
    }

    fn recompute_null(&mut self, enc: &Encoded) -> OperativeNull {
        self.null_recompute_count += 1;
        self.last_recompute_ncols = Some(enc.ncols);
        println!(
            "Null recompute #{} on pruned alignment ({} columns); resuming pruning.",
            self.null_recompute_count, enc.ncols
        );
        self.compute_null_state(enc)
    }

    fn stats(&self, raw_counts: &[i64], nrows: usize) -> Stats {
        let counts = chi2::fudge_counts(raw_counts.to_vec());
        let expected = chi2::expected_observed(&counts, self.n_chars, nrows);
        let per_row = chi2::row_chi2(&expected, &counts, self.n_chars, nrows);
        let sum = crate::npy::sum(&per_row);
        Stats { counts, expected, per_row, sum }
    }
}

fn reported_stats(per_row: &[f64], sums: &RobustStats, sums_raw: &[f64], threshold: f64) -> ReportedStats {
    let total = crate::npy::sum(per_row);
    ReportedStats {
        z: sums.zscore(total),
        p: permutation::empirical_p(total, sums_raw),
        biased: per_row.iter().filter(|&&v| v > threshold).count(),
    }
}

/// Undo the most recent batch (`_do_touchdown_rollback`). Returns
/// (encoded alignment, string-column map, original_indices, n_recent).
fn touchdown_rollback(
    top_n_indices: &[usize],
    prune_dict: &mut IndexMap<usize, PruneRow>,
    original: &Encoded,
) -> (Encoded, Vec<usize>, Vec<usize>, usize) {
    let n_recent = top_n_indices.len();
    let keys: Vec<usize> = prune_dict.keys().copied().collect();
    for k in &keys[keys.len() - n_recent..] {
        prune_dict.shift_remove(k);
    }
    let mut removed: Vec<usize> = prune_dict.keys().copied().collect();
    removed.sort_unstable();
    // Re-encoded from the string alignment: C-ordered again.
    let enc = original.delete_columns_raw(&removed);
    let mut is_removed = vec![false; original.ncols];
    for &c in &removed {
        is_removed[c] = true;
    }
    let indices: Vec<usize> = (0..original.ncols).filter(|&i| !is_removed[i]).collect();
    (enc, indices.clone(), indices, n_recent)
}

/// `AlignmentPruner.recursive_prune` on an encoded alignment.
fn recursive_prune(p: &mut Pruner, original: &Encoded, initial_null: OperativeNull) -> Result<PruneResult> {
    let opts = p.opts;
    let t_n = original.nrows;
    let alignment_size = original.ncols;
    let mut prune_dict: IndexMap<usize, PruneRow> = IndexMap::new();
    let mut iteration = 0usize;
    let mut removed_columns_count: i64 = 0;
    let mut original_indices: Vec<usize> = (0..alignment_size).collect();

    // Fixed references for the reported before/after summary.
    let original_sums_raw = initial_null.state.null.sums.clone();
    let original_sums = RobustStats::of(&original_sums_raw);
    let original_threshold = initial_null.state.biased_threshold;

    let mut null = initial_null;
    let mut sums_stats = original_sums;

    let mut top_n_indices: Vec<usize> = Vec::new();
    let mut initial_global_chi2: Option<f64> = None;
    let mut chi2_differences: Vec<f64> = Vec::new();
    let mut null_pvals_for_batch: HashMap<usize, f64> = HashMap::new();

    let mut in_touchdown = false;
    p.null_recompute_count = 0;
    p.last_recompute_ncols = None;

    let mut enc = original.clone();
    // Original column index of each column of the current (string) alignment.
    let mut cols: Vec<usize> = (0..alignment_size).collect();
    let mut count_rows_raw = chi2::count_rows_raw(&enc, p.n_chars);

    let mut stop_reason: Option<String> = None;
    let mut before_real: Option<Vec<f64>> = None;
    let mut before_stats: Option<ReportedStats> = None;
    let mut after_stats: Option<ReportedStats> = None;
    let mut per_row_chi2: Option<Vec<f64>> = None;

    while removed_columns_count + (top_n_indices.len() as i64) < opts.max_residue_pruned {
        let stats = p.stats(&count_rows_raw, t_n);
        let alignment_empirical_p = permutation::empirical_p(stats.sum, &null.state.null.sums);
        let significant_count =
            stats.per_row.iter().filter(|&&v| v > null.state.biased_threshold).count();
        let alignment_z = sums_stats.zscore(stats.sum);

        if iteration == 0 {
            before_real = Some(stats.per_row.clone());
            before_stats = Some(reported_stats(&stats.per_row, &original_sums, &original_sums_raw, original_threshold));
        } else {
            for &col in &top_n_indices {
                let pruned_col = original_indices[col];
                prune_dict.insert(
                    pruned_col,
                    PruneRow {
                        iteration,
                        original_index: pruned_col,
                        global_chi2: stats.sum,
                        initial_score: initial_global_chi2.expect("score computed before pruning"),
                        difference: chi2_differences[col],
                        significant_taxa: significant_count,
                        null_pvalue: null_pvals_for_batch.get(&col).copied(),
                    },
                );
                removed_columns_count += 1;
            }
            let mut drop = vec![false; original_indices.len()];
            for &c in &top_n_indices {
                drop[c] = true;
            }
            original_indices = original_indices
                .iter()
                .enumerate()
                .filter(|(j, _)| !drop[*j])
                .map(|(_, &i)| i)
                .collect();
        }
        per_row_chi2 = Some(stats.per_row.clone());

        let threshold = if opts.delta_null { 0.95 } else { 0.05 };
        let should_stop = alignment_empirical_p > threshold;

        println!(
            "Columns removed: {}, {}% | Biased taxa permutation: {} | Alignment Z-score: {} | Alignment p-value: {}",
            removed_columns_count,
            fixed(removed_columns_count as f64 / alignment_size as f64 * 100.0, 2),
            significant_count,
            fixed(alignment_z, 2),
            fixed(alignment_empirical_p, 2)
        );

        let target_top_n = (p.initial_top_n / 10).max(1);
        if should_stop {
            let mut rolled = false;
            if !in_touchdown && p.top_n > target_top_n && !top_n_indices.is_empty() {
                let (e, c, oi, n_recent) = touchdown_rollback(&top_n_indices, &mut prune_dict, original);
                enc = e;
                cols = c;
                original_indices = oi;
                count_rows_raw = chi2::count_rows_raw(&enc, p.n_chars);
                p.top_n = target_top_n;
                removed_columns_count -= n_recent as i64;
                in_touchdown = true;
                top_n_indices.clear();
                rolled = true;
                println!(
                    "Touchdown: alignment-level overshoot (reason='alignment p-value', alignment_p={}). \
                     Rolled back {n_recent} columns, reduced top_n to {target_top_n}, resuming pruning.",
                    fixed(alignment_empirical_p, 3)
                );
            }
            if p.recompute_due(enc.ncols, iteration) {
                null = p.recompute_null(&enc);
                sums_stats = RobustStats::of(&null.state.null.sums);
                top_n_indices.clear();
                continue;
            }
            if rolled {
                continue;
            }
            stop_reason = Some("alignment p-value".to_string());
            after_stats = Some(reported_stats(&stats.per_row, &original_sums, &original_sums_raw, original_threshold));
            break;
        }

        let (initial, deltas) = delta::score_and_deltas(
            p.alg,
            &stats.counts,
            &stats.expected,
            &enc,
            p.n_chars,
            Some(&null.z),
            opts.num_workers_chisq,
            true,
        );
        initial_global_chi2 = Some(initial);
        chi2_differences = deltas;

        // Stable sort by descending delta (Python sorted(..., reverse=True)).
        let mut order: Vec<usize> = (0..chi2_differences.len()).collect();
        order.sort_by(|&a, &b| {
            chi2_differences[b].partial_cmp(&chi2_differences[a]).unwrap_or(std::cmp::Ordering::Equal)
        });
        order.truncate(p.top_n);
        let sorted_candidates = order;

        let justified: Vec<usize>;
        if let Some(null_max) = &null.max_deltas {
            let mut just = Vec::new();
            null_pvals_for_batch.clear();
            let mut last_p = f64::NAN;
            for &col in &sorted_candidates {
                let d = chi2_differences[col];
                last_p = null_max.iter().filter(|&&v| v >= d).count() as f64 / null_max.len() as f64;
                if last_p > 0.05 {
                    break;
                }
                just.push(col);
                null_pvals_for_batch.insert(col, last_p);
            }
            if just.is_empty() {
                if p.recompute_due(enc.ncols, iteration) {
                    null = p.recompute_null(&enc);
                    sums_stats = RobustStats::of(&null.state.null.sums);
                    top_n_indices.clear();
                    continue;
                }
                let Some(&first) = sorted_candidates.first() else {
                    bail!("no candidate columns left to test");
                };
                let max_delta = chi2_differences[first];
                stop_reason = Some(format!(
                    "delta not significant (p={}, max_delta={})",
                    fixed(last_p, 3),
                    fixed(max_delta, 6)
                ));
                after_stats = Some(reported_stats(&stats.per_row, &original_sums, &original_sums_raw, original_threshold));
                break;
            }
            if just.len() < sorted_candidates.len() {
                let fail_delta = chi2_differences[sorted_candidates[just.len()]];
                let mut rolled = false;
                if !in_touchdown && p.top_n > target_top_n && !top_n_indices.is_empty() {
                    let (e, c, oi, n_recent) = touchdown_rollback(&top_n_indices, &mut prune_dict, original);
                    enc = e;
                    cols = c;
                    original_indices = oi;
                    count_rows_raw = chi2::count_rows_raw(&enc, p.n_chars);
                    p.top_n = target_top_n;
                    removed_columns_count -= n_recent as i64;
                    in_touchdown = true;
                    top_n_indices.clear();
                    rolled = true;
                    println!(
                        "Touchdown: delta-null partial batch (rank {} fails at p={}, delta={}). \
                         Rolled back {n_recent} columns, reduced top_n to {target_top_n}, resuming pruning.",
                        just.len() + 1,
                        fixed(last_p, 3),
                        fixed(fail_delta, 6)
                    );
                }
                if p.recompute_due(enc.ncols, iteration) {
                    null = p.recompute_null(&enc);
                    sums_stats = RobustStats::of(&null.state.null.sums);
                    top_n_indices.clear();
                    continue;
                }
                if rolled {
                    continue;
                }
                println!(
                    "  Partial batch: removing {}/{} columns (rank {} fails at p={}, delta={})",
                    just.len(),
                    sorted_candidates.len(),
                    just.len() + 1,
                    fixed(last_p, 3),
                    fixed(fail_delta, 6)
                );
            }
            justified = just;
        } else {
            justified = sorted_candidates;
        }

        top_n_indices = justified;
        for &c in &top_n_indices {
            for (t, &code) in enc.column(c).iter().enumerate() {
                if code >= 0 {
                    count_rows_raw[code as usize * t_n + t] -= 1;
                }
            }
        }
        enc = enc.delete_columns(&top_n_indices);
        let mut drop = vec![false; cols.len()];
        for &c in &top_n_indices {
            drop[c] = true;
        }
        cols = cols.iter().enumerate().filter(|(j, _)| !drop[*j]).map(|(_, &c)| c).collect();
        iteration += 1;
    }

    let cap_reason = format!("reached --max_residue_pruned cap ({})", opts.max_residue_pruned);
    if stop_reason.is_none() && !top_n_indices.is_empty() {
        let final_stats = p.stats(&count_rows_raw, t_n);
        let final_biased =
            final_stats.per_row.iter().filter(|&&v| v > null.state.biased_threshold).count();
        for &col in &top_n_indices {
            let pruned_col = original_indices[col];
            prune_dict.insert(
                pruned_col,
                PruneRow {
                    iteration,
                    original_index: pruned_col,
                    global_chi2: final_stats.sum,
                    initial_score: initial_global_chi2.expect("score computed before pruning"),
                    difference: chi2_differences[col],
                    significant_taxa: final_biased,
                    null_pvalue: null_pvals_for_batch.get(&col).copied(),
                },
            );
        }
        stop_reason = Some(cap_reason);
        after_stats = Some(reported_stats(&final_stats.per_row, &original_sums, &original_sums_raw, original_threshold));
        per_row_chi2 = Some(final_stats.per_row);
    } else if stop_reason.is_none() {
        stop_reason = Some(cap_reason);
        if after_stats.is_none() {
            after_stats = before_stats;
        }
    }

    let Some(after_real) = per_row_chi2 else {
        // Python crashes here (UnboundLocalError) when the loop never runs.
        bail!("--max_residue_pruned must be at least 1 (got {})", opts.max_residue_pruned);
    };
    p.top_n = p.initial_top_n;
    Ok(PruneResult {
        kept_columns: cols,
        prune_dict,
        before_real: before_real.expect("first iteration ran"),
        after_real,
        stop_reason: stop_reason.unwrap(),
        null_recomputes: p.null_recompute_count,
        null_max_deltas: null.max_deltas,
        before_stats: before_stats.expect("first iteration ran"),
        after_stats: after_stats.expect("after stats set"),
    })
}

/// `AlignmentPruner.run`: read, prune, write the four output files.
pub fn run(opts: &PruneOptions) -> Result<()> {
    let start = Instant::now();
    let Some(alg) = Algorithm::parse(&opts.pruning_algorithm) else {
        bail!("Unknown pruning algorithm: {}", opts.pruning_algorithm);
    };
    if opts.top_n == 0 {
        bail!("--top_n must be at least 1");
    }
    if opts.permutations == 0 {
        bail!("--permutations must be at least 1");
    }
    if opts.max_residue_pruned < 1 {
        bail!("--max_residue_pruned must be at least 1 (got {})", opts.max_residue_pruned);
    }
    crate::blas::init_single_threaded();
    println!("Starting pruning workflow");
    println!("Chi² calculator using {} worker(s)", opts.num_workers_chisq);
    println!("Permutation test using {} worker(s)", opts.num_workers_permute);
    println!("Pruning algorithm: {}", opts.pruning_algorithm);
    println!("Initial top_n: {}", opts.top_n);

    let aln: Alignment = alignment::read_alignment(Path::new(&opts.file), &opts.format)?;
    let seq_type = alignment::detect_sequence_type(&aln)?;
    println!("{}", seq_type.message());
    let char_set = seq_type.char_set();
    let enc = Encoded::from_rows(&aln.data, aln.nrows(), aln.ncols, char_set);

    let mut pruner = Pruner {
        opts,
        alg,
        n_chars: char_set.len(),
        top_n: opts.top_n,
        initial_top_n: opts.top_n,
        null_recompute_count: 0,
        last_recompute_ncols: None,
    };
    println!("Running {} permutations.", opts.permutations);
    println!("Using {} worker(s) for permutation", opts.num_workers_permute);
    let null = pruner.compute_null_state(&enc);
    let original_null = null.state.null.clone();
    let result = recursive_prune(&mut pruner, &enc, null)?;

    let files = output::write_all(opts, alg, &aln, &original_null, &result)?;

    let n_taxa = aln.nrows();
    let n_input = aln.ncols;
    let n_removed = result.prune_dict.len();
    let pct = if n_input > 0 { n_removed as f64 / n_input as f64 * 100.0 } else { 0.0 };
    let (b, a) = (result.before_stats, result.after_stats);
    println!(
        "Pruned {n_input} → {} columns (-{n_removed}, {}%) | Z {} → {} | p {} → {} | biased {}/{n_taxa} → {}/{n_taxa} | stop: {} | null recomputes: {}",
        n_input - n_removed,
        fixed(pct, 2),
        fixed(b.z, 2),
        fixed(a.z, 2),
        fixed(b.p, 2),
        fixed(a.p, 2),
        b.biased,
        a.biased,
        result.stop_reason,
        result.null_recomputes
    );
    println!("Pruned alignment saved to {}", files.fasta);
    println!("Columns pruned in order saved to: {}", files.tsv);
    println!("Taxa before/after Z and p-values printed to {}", files.scores);
    println!("Chiscore values from permutation testing printed to {}", files.json);
    println!("Execution time: {} seconds", fixed(start.elapsed().as_secs_f64(), 2));
    let _ = SequenceType::Dna;
    Ok(())
}
