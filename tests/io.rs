//! Alignment reading parity: every input under tests/io/inputs and tests/data
//! must read, detect and re-serialise exactly as Python WitChi does
//! (expected outputs from `scripts/make_io_fixtures.py`).

use std::path::Path;
use witchi::alignment::{SequenceType, detect_sequence_type, read_alignment, render_fasta};

fn describe(path: &Path) -> String {
    let fmt = if path.extension().is_some_and(|e| e == "nex") { "nexus" } else { "fasta" };
    let aln = match read_alignment(path, fmt) {
        Ok(a) => a,
        Err(_) => return "ERROR".to_string(),
    };
    let kind = match detect_sequence_type(&aln) {
        Ok(SequenceType::Dna) => "dna".to_string(),
        Ok(SequenceType::Protein) => "protein".to_string(),
        Ok(SequenceType::UnexpectedProtein(c)) => {
            format!("unexpected:{}", String::from_utf8(c).unwrap())
        }
        Err(_) => "detect-error".to_string(),
    };
    let mut out = format!("OK {kind}\n").into_bytes();
    let rows: Vec<Vec<u8>> = (0..aln.nrows()).map(|t| aln.row(t)).collect();
    render_fasta(&aln.records, rows.iter().map(Vec::as_slice), &mut out);
    String::from_utf8(out).unwrap()
}

#[test]
fn readers_match_python() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let rust_only: Vec<String> = std::fs::read_to_string(root.join("io/inputs/RUST_ONLY_ERRORS"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    let mut failures = Vec::new();
    let mut n = 0;
    for dir in ["io/inputs", "data"] {
        for entry in std::fs::read_dir(root.join(dir)).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_string_lossy().to_string();
            if name == "RUST_ONLY_ERRORS" {
                continue;
            }
            let expected =
                std::fs::read_to_string(root.join("io/expected").join(format!("{name}.txt")))
                    .unwrap_or_else(|_| panic!("no expected output for {name}"));
            let got = describe(&path);
            let ok = if rust_only.contains(&name) || expected.starts_with("ERROR") {
                got == "ERROR"
            } else {
                got == expected
            };
            if !ok {
                let first = |s: &str| s.lines().take(3).collect::<Vec<_>>().join(" | ");
                failures.push(format!("{name}: expected {:?}, got {:?}", first(&expected), first(&got)));
            }
            n += 1;
        }
    }
    assert!(failures.is_empty(), "{} of {n} inputs differ:\n{}", failures.len(), failures.join("\n"));
}
