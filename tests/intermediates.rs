//! Stage-by-stage bit-exact comparison with Python WitChi internals
//! (dumps from `scripts/dump_intermediates.py`, cases in tests/intermediates.tsv).

use std::collections::HashMap;
use std::path::Path;

use witchi::alignment::{detect_sequence_type, read_alignment};
use witchi::chi2::{self, Encoded};
use witchi::delta::{self, Algorithm, ZParams};
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

#[test]
fn deltas_match_python() {
    witchi::blas::init_single_threaded();
    let mut failures = Vec::new();
    for case in cases() {
        let d = load_dump(&case.name);
        let s = setup(&case);
        let alg = Algorithm::parse(&case.algorithm).unwrap();
        let st = NullState::new(NullMatrix::compute(&s.enc, s.n_chars, case.permutations, 4));
        let z = ZParams::from_null(&st);
        let counts = chi2::fudge_counts(chi2::count_rows_raw(&s.enc, s.n_chars));
        let e = chi2::expected_observed(&counts, s.n_chars, s.enc.nrows);
        let (initial, deltas) =
            delta::score_and_deltas(alg, &counts, &e, &s.enc, s.n_chars, Some(&z), 1, true);
        check(&case.name, "initial_score", &[initial], &d["initial_score"], &mut failures);
        check(&case.name, "deltas", &deltas, &d["deltas"], &mut failures);
        let nmd = delta::null_max_deltas(alg, &s.enc, s.n_chars, case.permutations, Some(&z), 4);
        check(&case.name, "null_max_deltas", &nmd, &d["null_max_deltas"], &mut failures);
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// Ad-hoc comparison against an external dump:
/// WITCHI_DUMP=<dump.txt> WITCHI_DUMP_INPUT=<alignment> WITCHI_DUMP_ARGS="--pruning_algorithm squared --permutations 20"
#[test]
#[ignore = "manual: needs WITCHI_DUMP / WITCHI_DUMP_INPUT"]
fn external_dump_matches_python() {
    witchi::blas::init_single_threaded();
    let dump = std::env::var("WITCHI_DUMP").expect("WITCHI_DUMP");
    let input = std::env::var("WITCHI_DUMP_INPUT").expect("WITCHI_DUMP_INPUT");
    let args = std::env::var("WITCHI_DUMP_ARGS").unwrap_or_default();
    let a: Vec<&str> = args.split_whitespace().collect();
    let get = |k: &str, d: &str| a.iter().position(|x| *x == k).map(|i| a[i + 1].to_string()).unwrap_or(d.to_string());
    let alg = Algorithm::parse(&get("--pruning_algorithm", "wasserstein")).unwrap();
    let perms: usize = get("--permutations", "100").parse().unwrap();
    let d: HashMap<String, Vec<f64>> = std::fs::read_to_string(&dump)
        .unwrap()
        .lines()
        .map(|l| {
            let mut it = l.split(' ');
            let key = it.next().unwrap().to_string();
            it.next();
            (key, it.map(|h| f64::from_bits(u64::from_str_radix(h, 16).unwrap())).collect())
        })
        .collect();
    let aln = read_alignment(Path::new(&input), &get("--format", "fasta")).unwrap();
    let cs = detect_sequence_type(&aln).unwrap().char_set();
    let enc = Encoded::from_rows(&aln.data, aln.nrows(), aln.ncols, cs);
    let st = NullState::new(NullMatrix::compute(&enc, cs.len(), perms, 8));
    let z = ZParams::from_null(&st);
    let counts = chi2::fudge_counts(chi2::count_rows_raw(&enc, cs.len()));
    let e = chi2::expected_observed(&counts, cs.len(), enc.nrows);
    let (initial, deltas) = delta::score_and_deltas(alg, &counts, &e, &enc, cs.len(), Some(&z), 1, true);
    let mut failures = Vec::new();
    check("ext", "null_pooled", &st.null.pooled, &d["null_pooled"], &mut failures);
    check("ext", "initial_score", &[initial], &d["initial_score"], &mut failures);
    check("ext", "deltas", &deltas, &d["deltas"], &mut failures);
    let bad: Vec<usize> = (0..deltas.len()).filter(|&i| deltas[i].to_bits() != d["deltas"][i].to_bits()).collect();
    eprintln!("mismatching delta columns (first 40): {:?}", &bad[..bad.len().min(40)]);
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}
