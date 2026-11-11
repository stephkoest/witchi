//! Stage-by-stage bit-exact comparison with Python WitChi internals
//! (dumps from `scripts/dump_intermediates.py`, cases in tests/intermediates.tsv).

use std::collections::HashMap;
use std::path::Path;

use witchi::alignment::{detect_sequence_type, read_alignment};
use witchi::chi2::{self, Encoded};
use witchi::permutation::{NullMatrix, NullState};

struct Case {
    name: String,
    input: String,
    format: String,
    algorithm: String,
    permutations: usize,
}

fn cases() -> Vec<Case> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let text = std::fs::read_to_string(root.join("intermediates.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let args: Vec<&str> = f[2].split_whitespace().collect();
            let get = |k: &str, d: &str| {
                args.iter().position(|a| *a == k).map(|i| args[i + 1].to_string()).unwrap_or(d.to_string())
            };
            Case {
                name: f[0].to_string(),
                input: f[1].to_string(),
                format: get("--format", "fasta"),
                algorithm: get("--pruning_algorithm", "wasserstein"),
                permutations: get("--permutations", "100").parse().unwrap(),
            }
        })
        .collect()
}

fn load_dump(name: &str) -> HashMap<String, Vec<f64>> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/intermediates/{name}.txt"));
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| {
            let mut it = l.split(' ');
            let key = it.next().unwrap().to_string();
            let n: usize = it.next().unwrap().parse().unwrap();
            let v: Vec<f64> = it.map(|h| f64::from_bits(u64::from_str_radix(h, 16).unwrap())).collect();
            assert_eq!(v.len(), n);
            (key, v)
        })
        .collect()
}

/// Compare bit patterns; report the first few mismatches.
fn check(case: &str, what: &str, got: &[f64], expected: &[f64], failures: &mut Vec<String>) {
    if got.len() != expected.len() {
        failures.push(format!("{case}/{what}: length {} vs expected {}", got.len(), expected.len()));
        return;
    }
    let bad: Vec<usize> =
        (0..got.len()).filter(|&i| got[i].to_bits() != expected[i].to_bits()).collect();
    if !bad.is_empty() {
        let i = bad[0];
        failures.push(format!(
            "{case}/{what}: {} of {} differ; first at [{i}]: got {:e} expected {:e}",
            bad.len(),
            got.len(),
            got[i],
            expected[i]
        ));
    }
}

pub struct Setup {
    pub enc: Encoded,
    pub n_chars: usize,
}

fn setup(case: &Case) -> Setup {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(&case.input);
    let aln = read_alignment(&path, &case.format).unwrap();
    let st = detect_sequence_type(&aln).unwrap();
    let cs = st.char_set();
    Setup { enc: Encoded::from_rows(&aln.data, aln.nrows(), aln.ncols, cs), n_chars: cs.len() }
}

#[test]
fn null_state_matches_python() {
    let mut failures = Vec::new();
    for case in cases() {
        let d = load_dump(&case.name);
        let s = setup(&case);
        let null = NullMatrix::compute(&s.enc, s.n_chars, case.permutations, 4);
        let st = NullState::new(null);
        check(&case.name, "null_pooled", &st.null.pooled, &d["null_pooled"], &mut failures);
        check(&case.name, "null_sums", &st.null.sums, &d["null_sums"], &mut failures);
        check(&case.name, "z_mean", &[st.z_mean], &d["z_mean"], &mut failures);
        check(&case.name, "z_scale", &[st.z_scale], &d["z_scale"], &mut failures);
        check(&case.name, "z_quantiles", &st.z_quantiles, &d["z_quantiles"], &mut failures);
        check(&case.name, "biased_threshold", &[st.biased_threshold], &d["biased_threshold"], &mut failures);

        let real = chi2::per_row_chi2_of(&s.enc, s.n_chars);
        check(&case.name, "real_per_row_chi2", &real, &d["real_per_row_chi2"], &mut failures);
        check(&case.name, "real_sum", &[witchi::npy::sum(&real)], &d["real_sum"], &mut failures);
        let _ = &case.algorithm;
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}
