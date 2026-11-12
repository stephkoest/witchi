//! Golden tests: run the Rust `witchi prune` on each scenario in a manifest and
//! byte-compare its output files against the Python reference outputs produced
//! by `scripts/make_golden.sh`.
//!
//! * Small tier: `tests/cases.tsv`, inputs `tests/data/`, refs `tests/golden/<case>/`.
//! * Large tier: `tests/cases_large.tsv`, inputs and refs under `tests/golden_large/`
//!   (not committed; generate with `scripts/make_golden.sh --large`).
//!
//! Both tests are `#[ignore]` until `prune` is implemented (PLAN.md phase 5):
//!     cargo test --release --test golden -- --ignored
//! Set `WITCHI_GOLDEN_FILTER=<substring>` to run a subset of cases, and
//! `WITCHI_GOLDEN_BIN=<path>` to test another binary (e.g. Python's `witchi`).
//! `stdout.txt` in a reference dir is informational and never compared.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

struct Case {
    name: String,
    input: String,
    args: Vec<String>,
}

fn read_manifest(path: &Path) -> Vec<Case> {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read manifest {}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let mut parts = l.splitn(3, '\t');
            let name = parts.next().unwrap().to_string();
            let input = parts
                .next()
                .unwrap_or_else(|| panic!("manifest line without input: {l:?}"))
                .to_string();
            let args = parts
                .next()
                .unwrap_or("")
                .split_whitespace()
                .map(str::to_string)
                .collect();
            Case { name, input, args }
        })
        .collect()
}

/// Describe the first difference between two byte strings by line.
fn first_difference(expected: &[u8], actual: &[u8]) -> String {
    let exp_lines: Vec<&[u8]> = expected.split(|&b| b == b'\n').collect();
    let act_lines: Vec<&[u8]> = actual.split(|&b| b == b'\n').collect();
    let show = |l: Option<&&[u8]>| match l {
        None => "<missing>".to_string(),
        Some(l) => {
            let s = String::from_utf8_lossy(l);
            let s = s.escape_debug().to_string();
            if s.len() > 160 { format!("{}…", &s[..160]) } else { s }
        }
    };
    for i in 0..exp_lines.len().max(act_lines.len()) {
        if exp_lines.get(i) != act_lines.get(i) {
            return format!(
                "first difference at line {}:\n      expected: {}\n      actual:   {}",
                i + 1,
                show(exp_lines.get(i)),
                show(act_lines.get(i))
            );
        }
    }
    "files differ (no line difference found)".to_string()
}

/// Run one case; returns a list of failure descriptions (empty on success).
fn run_case(case: &Case, data_dir: &Path, golden_dir: &Path, work_root: &Path) -> Vec<String> {
    let ref_dir = golden_dir.join(&case.name);
    if !ref_dir.is_dir() {
        return vec![format!("missing reference dir {}", ref_dir.display())];
    }
    let work = work_root.join(&case.name);
    let _ = fs::remove_dir_all(&work);
    fs::create_dir_all(&work).unwrap();
    fs::copy(data_dir.join(&case.input), work.join(&case.input))
        .unwrap_or_else(|e| panic!("copy input {}: {e}", case.input));

    // WITCHI_GOLDEN_BIN lets the harness be pointed at another implementation,
    // e.g. the Python `witchi` itself to check that the references reproduce.
    let bin = std::env::var("WITCHI_GOLDEN_BIN")
        .unwrap_or_else(|_| env!("CARGO_BIN_EXE_witchi").to_string());
    let output = Command::new(bin)
        .current_dir(&work)
        .arg("prune")
        .arg("--file")
        .arg(&case.input)
        .args(&case.args)
        .output()
        .expect("failed to spawn witchi binary");
    fs::write(work.join("rust_stdout.txt"), &output.stdout).unwrap();
    fs::write(work.join("rust_stderr.txt"), &output.stderr).unwrap();
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return vec![format!(
            "witchi exited with {} — stderr: {}",
            output.status,
            stderr.trim()
        )];
    }

    let mut failures = Vec::new();
    let mut refs: Vec<PathBuf> = fs::read_dir(&ref_dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.file_name().unwrap() != "stdout.txt")
        .collect();
    refs.sort();
    for expected_path in refs {
        let file_name = expected_path.file_name().unwrap();
        let actual_path = work.join(file_name);
        let expected = fs::read(&expected_path).unwrap();
        match fs::read(&actual_path) {
            Err(_) => failures.push(format!("{}: not written", file_name.to_string_lossy())),
            Ok(actual) if actual != expected => failures.push(format!(
                "{}: {}",
                file_name.to_string_lossy(),
                first_difference(&expected, &actual)
            )),
            Ok(_) => {}
        }
    }
    failures
}

fn run_tier(manifest: &Path, data_dir: &Path, golden_dir: &Path, tier: &str) {
    let filter = std::env::var("WITCHI_GOLDEN_FILTER").unwrap_or_default();
    let work_root = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("golden-{tier}"));
    let cases: Vec<Case> = read_manifest(manifest)
        .into_iter()
        .filter(|c| c.name.contains(&filter))
        .collect();
    assert!(!cases.is_empty(), "no cases selected (filter {filter:?})");

    let mut failed = Vec::new();
    for case in &cases {
        let failures = run_case(case, data_dir, golden_dir, &work_root);
        if failures.is_empty() {
            eprintln!("ok    {}", case.name);
        } else {
            eprintln!("FAIL  {}", case.name);
            for f in &failures {
                eprintln!("    {f}");
            }
            failed.push(case.name.clone());
        }
    }
    eprintln!(
        "{}/{} cases byte-identical; work dirs in {}",
        cases.len() - failed.len(),
        cases.len(),
        work_root.display()
    );
    assert!(failed.is_empty(), "failed cases: {}", failed.join(", "));
}

#[test]
#[ignore = "enable once `prune` is implemented (PLAN.md phase 5)"]
fn golden_small() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    run_tier(&root.join("cases.tsv"), &root.join("data"), &root.join("golden"), "small");
}

#[test]
#[ignore = "large tier: generate refs with scripts/make_golden.sh --large first"]
fn golden_large() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let golden = root.join("golden_large");
    if !golden.is_dir() {
        eprintln!("skipping: {} does not exist", golden.display());
        return;
    }
    run_tier(&root.join("cases_large.tsv"), &golden.join("data"), &golden, "large");
}

#[test]
fn manifest_is_consistent() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let cases = read_manifest(&root.join("cases.tsv"));
    let mut names = std::collections::HashSet::new();
    for c in &cases {
        assert!(names.insert(&c.name), "duplicate case name {}", c.name);
        assert!(root.join("data").join(&c.input).is_file(), "missing input {}", c.input);
        assert!(root.join("golden").join(&c.name).is_dir(), "missing refs for {}", c.name);
    }
}
