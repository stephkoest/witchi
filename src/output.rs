//! The four `witchi prune` output files, byte-compatible with Python.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};

use crate::alignment::{Alignment, Record, render_fasta};
use crate::chi2::{Decoder, Encoded};
use crate::delta::Algorithm;
use crate::permutation::{self, NullMatrix, RobustStats};
use crate::pruner::{PruneOptions, PruneResult};
use crate::pyfmt::{Json, json_float_list, json_string, newline_indent, repr_f64, tsv_row};
use crate::WITCHI_VERSION;

pub struct OutputFiles {
    pub fasta: String,
    pub tsv: String,
    pub scores: String,
    pub json: String,
}

/// Python `os.path.splitext(path)[0]`.
pub fn splitext_root(path: &str) -> &str {
    let sep = path.rfind('/').map(|i| i as isize).unwrap_or(-1);
    let Some(dot) = path.rfind('.') else {
        return path;
    };
    if (dot as isize) > sep {
        let mut i = (sep + 1) as usize;
        while i < dot {
            if path.as_bytes()[i] != b'.' {
                return &path[..dot];
            }
            i += 1;
        }
    }
    path
}

fn write(path: &str, bytes: &[u8]) -> Result<()> {
    fs::write(path, bytes).with_context(|| format!("cannot write {path}"))
}

fn create(path: &str) -> Result<BufWriter<File>> {
    let f = File::create(path).with_context(|| format!("cannot write {path}"))?;
    Ok(BufWriter::with_capacity(1 << 22, f))
}

/// Taxa per FASTA rendering task: 64 consecutive taxa share each cache line of
/// a column of the column-major codes.
const FASTA_ROWS_PER_TASK: usize = 64;

/// The pruned alignment as Biopython FASTA, streamed: groups of taxa are
/// gathered from the column-major codes and rendered in parallel, then
/// written in order.
fn write_pruned_fasta(
    path: &str,
    records: &[Record],
    enc: &Encoded,
    decoder: &Decoder,
    kept: &[usize],
) -> Result<()> {
    let mut w = create(path)?;
    let t_n = enc.nrows;
    let starts: Vec<usize> = (0..t_n).step_by(FASTA_ROWS_PER_TASK).collect();
    let batch = 2 * rayon::current_num_threads().max(1);
    for group in starts.chunks(batch) {
        let rendered: Vec<Vec<u8>> = group
            .par_iter()
            .map(|&t0| {
                let t1 = (t0 + FASTA_ROWS_PER_TASK).min(t_n);
                let (rows, l) = (t1 - t0, kept.len());
                let mut tile = vec![0u8; rows * l];
                for (jj, &j) in kept.iter().enumerate() {
                    for (r, &code) in enc.column(j)[t0..t1].iter().enumerate() {
                        tile[r * l + jj] = decoder.byte(code);
                    }
                }
                let mut out = Vec::with_capacity(rows * (l + l / 60 + 64));
                render_fasta(&records[t0..t1], (0..rows).map(|r| &tile[r * l..(r + 1) * l]), &mut out);
                out
            })
            .collect();
        for chunk in rendered {
            w.write_all(&chunk).with_context(|| format!("cannot write {path}"))?;
        }
    }
    w.flush().with_context(|| format!("cannot write {path}"))
}

/// The score dictionary as `json.dumps(doc, indent=4)`, streamed member by
/// member; the (P, T) permuted null is formatted in parallel, row by row.
fn write_score_json(
    path: &str,
    alg: Algorithm,
    records: &[Record],
    original_null: &NullMatrix,
    res: &PruneResult,
) -> Result<()> {
    let mut w = create(path)?;
    let put = |w: &mut BufWriter<File>, s: &str| -> Result<()> {
        w.write_all(s.as_bytes()).with_context(|| format!("cannot write {path}"))
    };
    let mut first = true;
    let mut key = |out: &mut String, k: &str| {
        if !first {
            out.push(',');
        }
        first = false;
        newline_indent(out, 1);
        json_string(out, k);
        out.push_str(": ");
    };
    let mut head = String::from("{");
    let small = [
        ("schema_version", Json::Int(3)),
        ("witchi_version", Json::Str(WITCHI_VERSION.to_string())),
        ("algorithm", Json::Str(alg.name().to_string())),
        ("stop_reason", Json::Str(res.stop_reason.clone())),
        ("null_recomputes", Json::Int(res.null_recomputes as i64)),
        ("taxa", Json::List(records.iter().map(|r| Json::Str(r.id.clone())).collect())),
    ];
    for (k, v) in &small {
        key(&mut head, k);
        v.write_indented(&mut head, 1);
    }
    key(&mut head, "before_real");
    json_float_list(&mut head, &res.before_real, 1);
    key(&mut head, "after_real");
    json_float_list(&mut head, &res.after_real, 1);
    key(&mut head, "before_permuted");
    put(&mut w, &head)?;

    let t = original_null.ntaxa;
    let p = original_null.permutations;
    if p == 0 {
        put(&mut w, "[]")?;
    } else {
        put(&mut w, "[")?;
        let batch = 2 * rayon::current_num_threads().max(1);
        let rows: Vec<usize> = (0..p).collect();
        for group in rows.chunks(batch) {
            let rendered: Vec<String> = group
                .par_iter()
                .map(|&i| {
                    let mut s = String::new();
                    if i > 0 {
                        s.push(',');
                    }
                    newline_indent(&mut s, 2);
                    json_float_list(&mut s, &original_null.pooled[i * t..(i + 1) * t], 2);
                    s
                })
                .collect();
            for s in rendered {
                put(&mut w, &s)?;
            }
        }
        let mut close = String::new();
        newline_indent(&mut close, 1);
        close.push(']');
        put(&mut w, &close)?;
    }

    let mut tail = String::new();
    key(&mut tail, "null_max_deltas");
    match res.null_max_deltas.as_deref() {
        Some(md) => json_float_list(&mut tail, md, 1),
        None => tail.push_str("null"),
    }
    newline_indent(&mut tail, 0);
    tail.push('}');
    put(&mut w, &tail)?;
    w.flush().with_context(|| format!("cannot write {path}"))
}

pub fn write_all(
    opts: &PruneOptions,
    alg: Algorithm,
    aln: &Alignment,
    enc: &Encoded,
    decoder: &Decoder,
    original_null: &NullMatrix,
    res: &PruneResult,
) -> Result<OutputFiles> {
    let root = splitext_root(&opts.file);
    let stem = format!("{root}_{}_s{}_pruned", alg.name(), opts.top_n);
    let files = OutputFiles {
        fasta: format!("{stem}.fasta"),
        tsv: format!("{stem}.tsv"),
        scores: format!("{stem}_scores.tsv"),
        json: format!("{stem}_score_dict.json"),
    };

    // Pruned alignment.
    write_pruned_fasta(&files.fasta, &aln.records, enc, decoder, &res.kept_columns)?;

    // Prune trace.
    let mut tsv = String::new();
    tsv_row(
        &mut tsv,
        &[
            "Iteration", "Original Index", "Global ChiScore", "Initial Significant ChiScore",
            "ChiScore Difference", "Significant taxa", "Null-Pvalue", "Algorithm",
        ],
    );
    for r in res.prune_dict.values() {
        let fields = [
            r.iteration.to_string(),
            r.original_index.to_string(),
            repr_f64(r.global_chi2),
            repr_f64(r.initial_score),
            repr_f64(r.difference),
            r.significant_taxa.to_string(),
            r.null_pvalue.map(repr_f64).unwrap_or_default(),
            alg.name().to_string(),
        ];
        let refs: Vec<&str> = fields.iter().map(String::as_str).collect();
        tsv_row(&mut tsv, &refs);
    }
    write(&files.tsv, tsv.as_bytes())?;

    // Per-taxon before/after Z and p (against the original null).
    let pooled = &original_null.pooled;
    let sorted_null = permutation::SortedNull::new(pooled);
    let p_before = sorted_null.p_per_taxon(&res.before_real);
    let p_after = sorted_null.p_per_taxon(&res.after_real);
    let rs = RobustStats::of(pooled);
    // dict keyed by taxon id: duplicates keep the first position, last values.
    let mut by_name: IndexMap<&str, [f64; 5]> = IndexMap::new();
    for (i, rec) in aln.records.iter().enumerate() {
        let bz = rs.zscore(res.before_real[i]);
        let az = rs.zscore(res.after_real[i]);
        by_name.insert(rec.id.as_str(), [bz, p_before[i], az, p_after[i], az - bz]);
    }
    let mut entries: Vec<(&str, [f64; 5])> = by_name.into_iter().collect();
    // sorted(..., key=before_zscore, reverse=True) is stable.
    entries.sort_by(|a, b| b.1[0].partial_cmp(&a.1[0]).unwrap_or(std::cmp::Ordering::Equal));
    let mut scores = String::new();
    tsv_row(&mut scores, &["Row", "Before-Z", "Before-Pvalue", "After-Z", "After-Pvalue", "Delta-Z"]);
    for (name, v) in entries {
        let f: Vec<String> = v.iter().map(|&x| repr_f64(x)).collect();
        tsv_row(&mut scores, &[name, &f[0], &f[1], &f[2], &f[3], &f[4]]);
    }
    write(&files.scores, scores.as_bytes())?;

    // Score dictionary JSON.
    write_score_json(&files.json, alg, &aln.records, original_null, res)?;
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::splitext_root;

    #[test]
    fn splitext_matches_python() {
        assert_eq!(splitext_root("a/b.faa"), "a/b");
        assert_eq!(splitext_root("b.tar.gz"), "b.tar");
        assert_eq!(splitext_root(".bashrc"), ".bashrc");
        assert_eq!(splitext_root("dir.d/file"), "dir.d/file");
        assert_eq!(splitext_root("..x"), "..x");
        assert_eq!(splitext_root("a/..b.c"), "a/..b");
        assert_eq!(splitext_root("noext"), "noext");
    }
}
