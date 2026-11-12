//! The four `witchi prune` output files, byte-compatible with Python.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::fs;

use crate::alignment::{Alignment, render_fasta};
use crate::delta::Algorithm;
use crate::permutation::{self, NullMatrix, RobustStats};
use crate::pruner::{PruneOptions, PruneResult};
use crate::pyfmt::{Json, repr_f64, tsv_row};
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

pub fn write_all(
    opts: &PruneOptions,
    alg: Algorithm,
    aln: &Alignment,
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
    let mut fasta = Vec::with_capacity(aln.nrows() * (res.kept_columns.len() + 64));
    let mut row = Vec::with_capacity(res.kept_columns.len());
    let rows: Vec<Vec<u8>> = (0..aln.nrows())
        .map(|t| {
            let full = aln.row(t);
            row.clear();
            row.extend(res.kept_columns.iter().map(|&j| full[j]));
            row.clone()
        })
        .collect();
    render_fasta(&aln.records, rows.iter().map(Vec::as_slice), &mut fasta);
    write(&files.fasta, &fasta)?;

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
    let p_before = permutation::empirical_p_per_taxon(&res.before_real, pooled);
    let p_after = permutation::empirical_p_per_taxon(&res.after_real, pooled);
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
    let floats = |v: &[f64]| Json::List(v.iter().map(|&x| Json::Float(x)).collect());
    let t = original_null.ntaxa;
    let before_permuted = Json::List(
        (0..original_null.permutations).map(|i| floats(&pooled[i * t..(i + 1) * t])).collect(),
    );
    let doc = Json::Obj(vec![
        ("schema_version".into(), Json::Int(3)),
        ("witchi_version".into(), Json::Str(WITCHI_VERSION.to_string())),
        ("algorithm".into(), Json::Str(alg.name().to_string())),
        ("stop_reason".into(), Json::Str(res.stop_reason.clone())),
        ("null_recomputes".into(), Json::Int(res.null_recomputes as i64)),
        (
            "taxa".into(),
            Json::List(aln.records.iter().map(|r| Json::Str(r.id.clone())).collect()),
        ),
        ("before_real".into(), floats(&res.before_real)),
        ("after_real".into(), floats(&res.after_real)),
        ("before_permuted".into(), before_permuted),
        (
            "null_max_deltas".into(),
            res.null_max_deltas.as_deref().map(floats).unwrap_or(Json::Null),
        ),
    ]);
    write(&files.json, doc.dumps_indent4().as_bytes())?;
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
