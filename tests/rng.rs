//! numpy RNG parity (fixtures from `scripts/make_fixtures.py`).

use witchi::rng::{Pcg64, permute_columns};

const FIX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

#[test]
fn pcg64_raw_matches_numpy() {
    let text = std::fs::read_to_string(format!("{FIX}/rng_raw.txt")).unwrap();
    for line in text.lines() {
        let mut parts = line.split(' ');
        let seed: u64 = parts.next().unwrap().parse().unwrap();
        let mut rng = Pcg64::from_seed(seed);
        for (i, hex) in parts.enumerate() {
            let expected = u64::from_str_radix(hex, 16).unwrap();
            assert_eq!(rng.next_u64(), expected, "seed {seed} draw {i}");
        }
    }
}

#[test]
fn column_permutation_matches_numpy() {
    let text = std::fs::read_to_string(format!("{FIX}/rng_perm.txt")).unwrap();
    let mut lines = text.lines().peekable();
    let mut checked = 0;
    while let Some(header) = lines.next() {
        let dims: Vec<usize> = header.strip_prefix("shape ").unwrap().split(' ').map(|v| v.parse().unwrap()).collect();
        let (rows, cols) = (dims[0], dims[1]);
        let read_matrix = |lines: &mut std::iter::Peekable<std::str::Lines>| {
            let columns: Vec<Vec<u8>> = (0..cols).map(|_| lines.next().unwrap().as_bytes().to_vec()).collect();
            let mut m = vec![0u8; rows * cols];
            for (j, c) in columns.iter().enumerate() {
                for t in 0..rows {
                    m[t * cols + j] = c[t];
                }
            }
            m
        };
        let original = read_matrix(&mut lines);
        while lines.peek().is_some_and(|l| l.starts_with("seed ")) {
            let seed: u64 = lines.next().unwrap()[5..].parse().unwrap();
            let expected = read_matrix(&mut lines);
            assert_eq!(permute_columns(&original, rows, cols, seed), expected, "shape {rows}x{cols} seed {seed}");
            checked += 1;
        }
    }
    assert_eq!(checked, 15);
}
