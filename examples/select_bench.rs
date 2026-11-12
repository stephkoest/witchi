//! Micro-benchmark: full sort vs rank selection on T values (development aid).
use std::time::Instant;
use witchi::npy::{self, QuantilePlan, SelectScratch};

fn main() {
    let mut s = 0x9e37_79b9_7f4a_7c15u64;
    let mut unit = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let t_n = 10_000;
    let plan = QuantilePlan::new(t_n, &npy::interior_positions(200));
    let ranks = plan.ranks().to_vec();
    let dists: Vec<(&str, Vec<f64>)> = vec![
        ("uniform", (0..t_n).map(|_| unit()).collect()),
        ("lognormal", (0..t_n).map(|_| { let u: f64 = (0..12).map(|_| unit()).sum::<f64>() - 6.0; (1.0 * u).exp() * 50.0 }).collect()),
        ("lognormal s=2", (0..t_n).map(|_| { let u: f64 = (0..12).map(|_| unit()).sum::<f64>() - 6.0; (2.0 * u).exp() * 50.0 }).collect()),
        ("bulk+outliers", (0..t_n).map(|i| if i % 100 == 0 { 1e6 * unit() } else { 50.0 + 10.0 * unit() }).collect()),
    ];
    let reps = 300;
    for (name, v) in &dists {
        let mut work = v.clone();
        let t0 = Instant::now();
        for _ in 0..reps {
            work.copy_from_slice(v);
            work.sort_unstable_by(f64::total_cmp);
        }
        let sort = t0.elapsed().as_secs_f64() / reps as f64;
        let mut out = vec![0.0; ranks.len()];
        let mut sc = SelectScratch::default();
        let t0 = Instant::now();
        for _ in 0..reps {
            work.copy_from_slice(v);
            npy::select_ranks(&work, &ranks, &mut out, &mut sc);
        }
        let sel = t0.elapsed().as_secs_f64() / reps as f64;
        let t0 = Instant::now();
        for _ in 0..reps {
            work.copy_from_slice(v);
        }
        let copy = t0.elapsed().as_secs_f64() / reps as f64;
        println!("{name:14} sort {:7.1} us  select {:7.1} us  (copy {:5.1} us)", sort * 1e6, sel * 1e6, copy * 1e6);
    }
}
