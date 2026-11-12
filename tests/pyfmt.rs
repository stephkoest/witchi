//! Formatting parity against fixtures from `scripts/make_fixtures.py`.

use witchi::pyfmt::{Json, fixed, repr_f64, tsv_row};

const FIX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

#[test]
fn float_repr_and_fixed_match_python() {
    let text = std::fs::read_to_string(format!("{FIX}/pyfmt_floats.tsv")).unwrap();
    let mut n = 0;
    let mut failures = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let x = f64::from_bits(u64::from_str_radix(f[0], 16).unwrap());
        let got = [repr_f64(x), fixed(x, 2), fixed(x, 3), fixed(x, 4), fixed(x, 6)];
        for (g, e) in got.iter().zip(&f[1..]) {
            if g != e && failures.len() < 20 {
                failures.push(format!("{} ({x:e}): expected {e:?}, got {g:?}", f[0]));
            }
        }
        n += 1;
    }
    assert!(failures.is_empty(), "{} mismatches, e.g.:\n{}", failures.len(), failures.join("\n"));
    assert!(n > 10_000);
}

#[test]
fn json_matches_python() {
    let expected = std::fs::read_to_string(format!("{FIX}/pyfmt_json.txt")).unwrap();
    let s = |v: &str| Json::Str(v.to_string());
    let fl = |v: f64| Json::Float(v);
    let doc = Json::Obj(vec![
        ("schema_version".into(), Json::Int(3)),
        ("witchi_version".into(), s("0.2.0")),
        ("stop_reason".into(), s("delta not significant (p=0.060, max_delta=8.542055)")),
        ("none".into(), Json::Null),
        ("tricky".into(), s("q\"b\\s/\n\r\t\u{8}\u{c}\u{1}\u{7f} é ✓ 𝄞 end")),
        (
            "floats".into(),
            Json::List(
                [0.0, -0.0, 1.5, 1e-07, 1e16, f64::INFINITY, f64::NEG_INFINITY, f64::NAN]
                    .into_iter()
                    .map(fl)
                    .collect(),
            ),
        ),
        (
            "nested".into(),
            Json::List(vec![
                Json::List(vec![fl(1.0), fl(2.5)]),
                Json::List(vec![]),
                Json::List(vec![fl(3.25)]),
            ]),
        ),
        ("empty_obj".into(), Json::Obj(vec![])),
        (
            "obj".into(),
            Json::Obj(vec![
                ("k".into(), Json::List(vec![Json::Int(1), Json::Int(2)])),
                ("z".into(), Json::Int(0)),
            ]),
        ),
    ]);
    assert_eq!(doc.dumps_indent4(), expected);
}

#[test]
fn csv_matches_python() {
    let expected = std::fs::read_to_string(format!("{FIX}/pyfmt_csv.txt")).unwrap();
    let mut out = String::new();
    tsv_row(&mut out, &["Row", "a b", "tab\there", "quote\"inside", "nl\nx", "cr\rx", "", "plain"]);
    tsv_row(&mut out, &[""]);
    tsv_row(&mut out, &["x", ""]);
    tsv_row(&mut out, &["'single'", "comma,ok", "ünï"]);
    assert_eq!(out, expected);
}
