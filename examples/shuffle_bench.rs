//! Micro-benchmark of the permutation inner loops (development aid):
//! shuffling a column and adding its residue counts.
use std::time::Instant;
use witchi::chi2;
use witchi::rng::{Pcg64, U32Stream};

fn main() {
    for &t_n in &[10_000usize, 200_000] {
        let c_n = 20;
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let col: Vec<i8> = (0..t_n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                let r = (s >> 40) % 24;
                if r >= 20 { -1 } else { r as i8 }
            })
            .collect();
        let cols = (20_000_000 / t_n).max(20);
        let mut work = col.clone();
        let mut rng = Pcg64::from_seed(12345);
        let t0 = Instant::now();
        for _ in 0..cols {
            work.copy_from_slice(&col);
            rng.shuffle(&mut work);
        }
        let shuffle = t0.elapsed().as_secs_f64() / (cols * t_n) as f64;
        let mut stream = U32Stream::from_seed(12345);
        let t0 = Instant::now();
        for _ in 0..cols {
            work.copy_from_slice(&col);
            stream.shuffle(&mut work);
        }
        let fast = t0.elapsed().as_secs_f64() / (cols * t_n) as f64;
        println!("T {t_n:7}: stream shuffle {:.2} ns/cell", fast * 1e9);
        let mut counts = vec![0i64; c_n * t_n];
        let t0 = Instant::now();
        for _ in 0..cols {
            chi2::add_column_counts(&mut counts, &work, t_n);
        }
        let count = t0.elapsed().as_secs_f64() / (cols * t_n) as f64;
        println!("T {t_n:7}: shuffle {:.2} ns/cell, count {:.2} ns/cell (checksum {})",
            shuffle * 1e9, count * 1e9, counts.iter().sum::<i64>());
    }
}
