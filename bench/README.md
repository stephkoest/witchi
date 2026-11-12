# Thread-scaling benchmark of `witchi prune`

Machine: 2 × AMD EPYC 7402 (48 cores, 96 hardware threads, 2 NUMA nodes), otherwise idle.
Binary: the phase 8 port (before the optimizations in `OPTIMIZATION.md`).

- **Inputs:** 9 subsets of the GTDB r237 ar53 MSA: 100 / 1,000 / 10,000 taxa ×
  100 / 1,000 / 10,000 columns, with taxa and columns evenly spaced
  (`scripts/make_bench_inputs.py`, written to `bench/data/`, not committed).
- **Settings:** WitChi's defaults: wasserstein, 100 permutations, `top_n` 1,
  `--max_residue_pruned` 100, delta-null on.
- **Grid:** worker counts 1, 2, 4, 8, 16, 32, 48 and 96, raising only
  `--num_workers_permute`, only `--num_workers_chisq`, or both.
  - Runs are sequential.
  - Configurations under 20 s were run 3 times; tables use the median.
- **Files:**
  - `scaling.tsv`: every run (`scripts/bench_scaling.py`).
  - `scaling_tables.md`: the tables (`scripts/bench_report.py`).

## Findings

- **Which flag is raised doesn't matter.** Permute-only, chisq-only and both give the same
  time at the same worker count: the port runs one thread pool sized by the larger of the
  two. `--num_workers_chisq` only changes Python's column chunking (T·L ≥ 5 × 10⁶). In this
  grid that never changed the output.
- **Output doesn't depend on worker settings.** Each dataset produced a single output
  checksum across all of its runs (522 runs in total).
- **Scaling grows with problem size.** At 96 threads the speed-up is 31.7× on
  10,000 × 10,000 (1,020 s → 32 s). It is 24× on 1,000 × 10,000 and 12× on 10,000 × 1,000.
  Inputs of 100 × 100 run in about 0.1 s and don't benefit.
- **Hyperthreads help a little at the largest size.** 48 → 96 threads gives 36.5 s → 32.1 s.
- **Limits seen:**
  - **Per-taxon p-values** ran on one thread and grew with taxa²: about 4.2 s at 10,000
    taxa, whatever the column count. This capped 10,000 × 100 at 3.4×. Fixed by T1.
  - **Delta-null memory** held one permuted copy per thread. Peak memory on
    10,000 × 10,000 grew from 0.55 GB at 1 thread to 10.8 GB at 96. Fixed by M1.
  - **Too few column blocks at 1,000 taxa.** The block width is 250, so 1,000 columns is
    4 blocks, and prune iterations stop speeding up beyond 4 threads (1.04 s from w=4
    on). See T7.
