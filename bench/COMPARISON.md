# Before vs after optimization

This compares two builds of the Rust port:
- **Before:** `84413af`, the complete port before any optimization. Its tree is the one released as v0.3.0, apart from the version number.
- **After:** `611c6bc`, `rust-port` with PRs #10–#16 merged (T1–T7, M1–M5).

Both builds use the same OpenBLAS and ran on the same machine: 2 × AMD EPYC 7402, 48 cores and 96 hardware threads.

## How it was run

- **Script:** `scripts/bench_compare.py bench/compare_cases.tsv OUT before=… after=… --reps 2`.
- **Alternation:** the two binaries alternate within each repetition, so load changes affect both.
- **Repetitions:** cases under a minute ran twice, and the table gives the minimum. The single-thread cases and the 100,000- and 200,000-taxa cases ran once.
- **Measurement:** wall time and peak RSS come from GNU `time`. Stage times come from `WITCHI_TIMING=1`.
- **Checks:** every run records the MD5 of the four output files and of stdout, without its "Execution time" line.

The results are in `compare_summary.tsv` (per case) and `compare_runs.tsv` (per run).

The inputs:

| Input | Source |
|---|---|
| `gtdb_t*_c*` | the 9 scaling-benchmark subsets of the GTDB r237 ar53 MSA (`scripts/make_bench_inputs.py`) |
| Full MSA | `test_data/gtdb_r237_ar53_concatenated.masked.faa`, 12,944 × 13,540 |
| 20,000 × 2,000 | the new large golden input (`scripts/make_inputs.py --large`) |
| 51,776 × 2,000 | the full MSA stacked four times, 2,000 evenly spaced columns |
| 100,000 and 200,000 × 2,000 | synthetic taxa: copies of source taxa with 15 % of their columns resampled from the same column (`scripts/make_scale_inputs.py`) |

## Output

**Every run gave identical output**: the same checksums for all four files and stdout, before and after, on all 18 cases. The large golden tier adds a second check. At 20,000 taxa, above the 15,625 where the column-block width stops shrinking, the optimized binary matches Python v0.2.0 byte for byte.

## Wall time and peak memory

Defaults are wasserstein, 100 permutations, `top_n` 1, a cap of 100 columns and the delta-null on. "Large" settings are `top_n` 10, 20 permutations and a cap of 50. Thread count is set with `--num_workers_chisq` and `--num_workers_permute`.

| Input (taxa × columns) | Settings | Threads | Wall before | Wall after | Speed-up | Peak RSS before | Peak RSS after |
|---|---|---:|---:|---:|---:|---:|---:|
| 100 × 100 | defaults | 32 | 0.06 s | 0.06 s | 1.0× | 75 MB | 53 MB |
| 100 × 1,000 | defaults | 32 | 0.19 s | 0.16 s | 1.2× | 92 MB | 84 MB |
| 100 × 10,000 | defaults | 32 | 0.55 s | 0.50 s | 1.1× | 127 MB | 103 MB |
| 1,000 × 100 | defaults | 32 | 0.31 s | 0.16 s | 1.9× | 148 MB | 138 MB |
| 1,000 × 1,000 | defaults | 32 | 1.63 s | 0.92 s | 1.8× | 248 MB | 221 MB |
| 1,000 × 10,000 | defaults | 32 | 4.69 s | 3.41 s | 1.4× | 573 MB | 244 MB |
| 10,000 × 100 | defaults | 32 | 5.13 s | 0.95 s | 5.4× | 479 MB | 541 MB |
| 10,000 × 1,000 | defaults | 32 | 9.80 s | 3.50 s | 2.8× | 805 MB | 534 MB |
| 10,000 × 10,000 | defaults | 32 | 46.7 s | 20.4 s | 2.3× | 3,745 MB | 575 MB |
| 10,000 × 10,000 | defaults | 96 | 32.1 s | 14.0 s | 2.3× | 10,733 MB | 1,176 MB |
| 1,000 × 1,000 | defaults | 1 | 9.74 s | 6.30 s | 1.5× | 91 MB | 91 MB |
| 1,000 × 10,000 | defaults | 1 | 90.4 s | 61.7 s | 1.5× | 72 MB | 29 MB |
| 10,000 × 1,000 | defaults | 1 | 103.8 s | 51.3 s | 2.0× | 128 MB | 51 MB |
| Full MSA, 12,944 × 13,540 | large | 32 | 18.1 s | 8.5 s | 2.1× | 4,164 MB | 564 MB |
| 20,000 × 2,000 | large | 32 | 7.52 s | 2.47 s | 3.0× | 1,264 MB | 616 MB |
| 51,776 × 2,000 | large | 32 | 30.5 s | 6.5 s | 4.7× | 3,205 MB | 1,395 MB |
| 100,000 × 2,000 | large | 32 | 124.6 s | 12.0 s | 10.4× | 6,176 MB | 2,363 MB |
| 200,000 × 2,000 | large | 32 | 1,110 s | 23.5 s | 47× | 12,313 MB | 4,886 MB |

### Where the time went

| Input | Threads | Build | Read | Null matrix | Delta-null | Prune iterations | Write outputs |
|---|---:|---|---:|---:|---:|---:|---:|
| 10,000 × 10,000 | 32 | before | 0.35 s | 4.29 s | 20.5 s | 13.4 s | 4.53 s |
| | | after | 0.19 s | 3.05 s | 10.5 s | 5.55 s | 0.46 s |
| Full MSA | 32 | before | 0.58 s | 2.49 s | 10.4 s | 1.32 s | 2.25 s |
| | | after | 0.33 s | 1.74 s | 4.86 s | 0.56 s | 0.56 s |
| 200,000 × 2,000 | 32 | before | 1.49 s | 6.84 s | 39.0 s | 7.27 s | 1,052 s |
| | | after | 0.82 s | 3.45 s | 12.1 s | 2.15 s | 2.58 s |

- **Large taxon counts:** the speed-up grows with taxa because the old per-taxon p-values were quadratic in the number of taxa (T1). At 200,000 taxa they took 1,050 of the 1,110 seconds, even with only 20 permutations.
- **Other stages:** the delta-null runs 2–3× faster and prune iterations 2.4–3.4× faster.
- **Memory:** peak memory no longer grows with the number of threads multiplied by the alignment size (M1). On 10,000 × 10,000 at 96 threads it went from 10.7 GB to 1.2 GB.
- **Small inputs:** inputs of 100 taxa run in under a second either way.
- **Slight regression:** 10,000 × 100 at 32 threads uses 62 MB more after (541 MB against 479 MB). This is the per-thread scratch (T4, T6), which outweighs the tiny alignment.

## Profile near the target size

The target is about 200,000 taxa × 40,000 columns. The optimized binary was run on a synthetic 200,000 × 10,000 alignment (2 GB) with WitChi's defaults on 48 threads. The machine's load average reached 45 from other users during the run, so these times are upper bounds.

| Stage | Time | Share | Scales with |
|---|---:|---:|---|
| Read alignment | 5.3 s | 1 % | cells |
| Null matrix (100 permutations) | 40.8 s | 10 % | cells × permutations |
| Null state | 1.9 s | 0.5 % | permutations × taxa |
| Delta-null (100 permutations) | 187 s | 46 % | cells × permutations |
| Prune iterations (100 × 1.56 s) | 156 s | 38 % | cells × iterations |
| Write outputs | 3.6 s | 1 % | cells |
| **Total wall** | **407 s** | | peak RSS 12.0 GB |

**Extrapolated to 200,000 × 40,000 on 48 threads:**
- **Time:** every large stage is linear in the column count, so roughly 4 × 400 s, about 27 minutes:
  - delta-null about 12.5 minutes
  - prune iterations about 10 minutes
  - null matrix about 3 minutes
- **Memory:** about 18 GB. That is 8 GB for the alignment plus the roughly 10 GB of per-thread working memory measured here, which depends on the taxon and thread counts, not on the column count.

The old binary could not have run this size:
- Its per-taxon p-values alone would take about 1.5 hours with 100 permutations.
- Its delta-null kept one 8 GB permuted copy per thread, which is 384 GB at 48 threads.

**Where the remaining time is:** the delta-null and the prune iterations are nearly all of it. Both are the shuffle plus the delta kernel (BLAS product, per-row block, Wasserstein rank selection), which the byte-identical work has already optimized. Further large gains would have to change how much work is done, and that changes the output: fewer delta-null permutations, a larger `top_n` (fewer iterations), or approximate quantiles.
