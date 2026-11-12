# Optimization review: byte-identical changes for very large alignments

Target: about 200,000 sequences × 40,000 amino-acid columns (8 × 10⁹ cells), default
`prune` settings (wasserstein, 100 permutations, delta-null on).

Every change listed here keeps the output files byte-identical. Each item says why.
Changes that would alter output (different block widths, different BLAS shapes, f32,
approximate quantiles, fewer permutations) are out of scope.

**Basis.** This is a code review plus extrapolation from the stage timings in
`bench/scaling.tsv` (alignments up to 10,000 × 10,000). Nothing was profiled at the target
size, so the per-cell costs and totals below are estimates.

Measured single-thread costs used for the estimates (1,000 × 10,000, wasserstein):

| Stage | Cost |
|---|---|
| Shuffle one column and add its counts (null matrix) | about 10 ns per cell |
| Delta kernel (prune iteration) | about 33 ns per cell at 1,000 taxa |
| Delta-null (shuffle + copy + kernel) | about 46 ns per cell |
| Per-taxon p-values at output | 4.2 s at 10,000 taxa, growing with taxa² |

## Status

Implemented: **T1, M2, M1** (PR #10), then **T2**. With them, all parity suites pass:
- the 33 + 3 + 9 goldens
- the 15 intermediate cases
- the I/O, formatting and RNG fixtures
- the full-MSA run, byte-identical to Python

Unit tests check, bit for bit:
- the sorted-null count against the elementwise count
- the streamed delta-null against the materialized one
- rank selection against a full sort (adversarial value sets, forced bucketing at small sizes)
- the Wasserstein reducer against the old sort-based reducer

### T1, M2, M1

Measured on the AMD EPYC 7402 (same runs as `bench/scaling.tsv`, old → new):

| Input | Threads | Wall (s) | Peak RSS (MB) | Notes |
|---|---:|---|---|---|
| 10,000 × 100 | 32 | 5.25 → 2.02 | 478 → 406 | output step 4.2 s → 0.9 s (T1) |
| 10,000 × 1,000 | 96 | 8.66 → 5.27 | 2,193 → 1,164 | |
| 10,000 × 10,000 | 32 | 46.4 → 44.6 | 3,751 → 841 | delta-null 20.5 s → 24.2 s (second shuffle pass) |
| 10,000 × 10,000 | 96 | 32.2 → 30.3 | 10,814 → 1,479 | |
| Full MSA 12,944 × 13,540 | 32 | 16.4 → 16.7 | 4,269 → 1,050 | top_n 10, 20 permutations, cap 50 |

### T2

Measured the same way (wall time; before → after T2):

| Input | Threads | Wall (s) | Delta-null (s) | Prune iterations (s) |
|---|---:|---|---|---|
| 10,000 × 1,000 | 1 (20 perms, cap 20) | – | 12.0 → 7.8 | 7.98 → 3.97 |
| 10,000 × 1,000 | 32 | 6.87 → 5.14 | 2.68 → 1.91 | 2.32 → 1.31 |
| 10,000 × 10,000 | 32 | 44.6 → 29.4 | 24.2 → 15.7 | 13.5 → 6.8 |
| 10,000 × 10,000 | 96 | 30.3 → 21.6 | 16.2 → 10.9 | 9.0 → 5.5 |
| Full MSA 12,944 × 13,540 | 32 | 16.7 → 12.5 | 11.0 → 7.4 | 1.28 → 0.73 |

Single-threaded, the squared algorithm (same kernel, no sort) takes 6.2 s and 2.2 s for those
two stages on 10,000 × 1,000. So T2 removed about 70 % of Wasserstein's extra cost; most of
what is left is the strided per-column gather (T3).

Next by impact: T3 and T6 (kernel loop order and buffers), then the I/O memory items M3–M5.

## What breaks first at 200,000 × 40,000

One copy of the alignment is 8 GB.

| Problem | Current behaviour | Estimate at target size |
|---|---|---|
| Delta-null memory | One full permuted copy per worker thread | 8 GB × threads: 384 GB at 48, 768 GB at 96 (this machine has 502 GB) |
| Base memory | Row-major chars + column-major codes + a clone + a new copy per iteration | 24–32 GB steady |
| Load memory | File text + per-record vectors + row-major + column-major | about 30 GB peak |
| Per-taxon p-values | Compares every taxon against the whole pooled null | about 28 minutes on one thread |
| Column deletion | Copies the whole encoded alignment every iteration | 8 GB allocated and copied per iteration |
| Delta kernel | Full sort per column; cache-hostile loop order once blocks exceed the cache | most of the run time |

With all the memory items done, the footprint should be about 8 GB plus roughly 0.2 GB
per worker thread.

## Memory

### M1. Delta-null: stop materializing the permuted alignment (highest priority)

- **Where:** [src/delta.rs:389](src/delta.rs:389) `null_max_deltas`,
  [src/permutation.rs:27](src/permutation.rs:27) `permuted`.
- **Now:** each permutation builds a complete permuted copy (T × L bytes) and all
  in-flight permutations hold theirs at once. In the 10,000 × 10,000 pilot this was
  3.8 GB at 32 workers for a 100 MB alignment.
- **Change:** two streaming passes per permutation. Pass 1 shuffles column by column and
  accumulates the per-taxon counts (exactly what `permuted_row_chi2` does). Pass 2
  restarts the generator from the same seed and regenerates the permuted columns one
  kernel block at a time (16 columns = 3.2 MB at 200,000 taxa), feeding each block to
  `per_row_block`.
- **Identical because:** the same seed gives the same permuted columns; the kernel sees
  the same values in the same blocks.
- **Cost:** the shuffle runs twice. Today that adds about 20 % to the delta-null stage
  (10 ns of 46 ns per cell); less after T4. If that matters on small inputs, keep the
  current path when `T × L × threads` is under a memory budget.

### M2. Prune loop: delete columns by index, not by copying

- **Where:** [src/pruner.rs:231](src/pruner.rs:231) `original.clone()`,
  [src/pruner.rs:437](src/pruner.rs:437) `enc.delete_columns`,
  [src/chi2.rs:45](src/chi2.rs:45), and the rollback at
  [src/pruner.rs:182](src/pruner.rs:182).
- **Now:** the loop starts from a clone of the encoded alignment and builds a new
  compacted copy every iteration and on rollback.
- **Change:** make `Encoded` a view: shared, immutable code matrix plus a vector of
  active column indices (the pruner already maintains this as `cols`). Deleting columns
  edits the index vector; rollback rebuilds it. Keep the `fortran` flag logic unchanged.
- **Identical because:** `column(j)` returns the same bytes; nothing numeric changes.
- **Effect:** removes 8–16 GB and an 8 GB allocate-and-copy per iteration. At the target
  size that copy is a single-threaded step of a few seconds next to an iteration that
  should take about 10 s on 48 cores.

### M3. Keep one copy of the alignment, not two

- **Where:** `Alignment.data` (row-major chars) and `Encoded.codes` (column-major codes),
  [src/chi2.rs:27](src/chi2.rs:27).
- **Now:** both live for the whole run; the chars are only needed again to write the
  pruned FASTA.
- **Change:** make the codes lossless. Residues in the char set keep codes 0..C−1; every
  other byte gets its own negative code (input is ASCII and upper-cased, so there are at
  most 124 such bytes and an `i8` has 128 negative values). All kernels already test
  `code >= 0`. Drop the row-major chars after encoding and decode at output time.
- **Identical because:** kernels see the same sign and the same non-negative codes; the
  output decodes to the same characters.
- **Effect:** steady-state alignment memory 16 GB → 8 GB.

### M4. Loading: stream the file

- **Where:** [src/alignment.rs:72](src/alignment.rs:72) `read_text`,
  [src/fasta.rs](src/fasta.rs), [src/alignment.rs:43](src/alignment.rs:43) `from_records`.
- **Now:** the whole file is read into a string, each record gets its own growable
  vector (capacity can be up to twice the sequence), then everything is copied into the
  row-major matrix and again into the codes.
- **Change:** read through a buffered reader, validate UTF-8 and translate newlines per
  chunk, and append residues straight into one buffer. Then transpose and encode into
  the column-major matrix (see T5) and free the buffer. A two-pass read (count records,
  then fill) avoids even the transient second copy.
- **Identical because:** same parsing rules, same bytes.
- **Effect:** load peak about 30 GB → 16 GB (8 GB with the two-pass read).

### M5. Output: stream the FASTA and the JSON

- **Where:** [src/output.rs:62](src/output.rs:62) (all pruned rows, then the whole FASTA
  in one buffer) and [src/output.rs:127](src/output.rs:127) (a JSON tree with one node
  per float: 20 million nodes for the permuted-null matrix).
- **Change:** write the FASTA through a buffered writer, transposing a tile of rows at a
  time from the column-major codes. Write the JSON incrementally, formatting each
  permutation's row of floats in parallel and emitting the rows in order.
- **Identical because:** same bytes in the same order.
- **Effect:** removes about 16 GB of transient memory for the FASTA and about 1 GB for
  the JSON; float formatting (20 million `repr` calls) drops from tens of seconds on one
  thread to about a second.

## Time

### T1. Per-taxon p-values are quadratic in the number of taxa (highest priority)

- **Where:** [src/permutation.rs:85](src/permutation.rs:85) `empirical_p_per_taxon`,
  called twice from [src/output.rs:100](src/output.rs:100) and once from
  [src/testcmd.rs:90](src/testcmd.rs:90).
- **Now:** for each taxon it scans the whole pooled null (P × T values): P × T²
  comparisons per call. The benchmark shows it: the output stage takes 0.02 s at 100
  taxa, 0.14 s at 1,000 and 4.2 s at 10,000, whatever the column count.
- **Change:** sort the pooled null once (NaNs removed), then binary-search each taxon:
  the count of null values ≥ x is the number of non-NaN values minus the number below x.
- **Identical because:** it is the same integer count, then the same division.
- **Effect:** about 28 minutes → a few seconds at 200,000 taxa. Also removes most of
  the fixed 4 s at 10,000 taxa, which currently caps thread scaling there.

### T2. Wasserstein: select the needed order statistics, don't sort

- **Where:** [src/delta.rs:265](src/delta.rs:265), full sort at
  [src/delta.rs:275](src/delta.rs:275).
- **Now:** every column's T z-scores are fully sorted, but only 2 × 200 fixed ranks are
  read (`QuantilePlan` prev/next). The sort is most of the kernel.
- **Done** (`npy::select_ranks`, `QuantilePlan::ranks` / `apply_ranked`):
  - **Selection on raw values:** the reducer selects the plan's ranks on the raw per-row
    values and converts only those to z-scores. This needs a finite mean and a positive,
    finite scale; otherwise, or when a column holds NaN, it behaves as before.
  - **Bucketing by bit pattern:** values are mapped to their total-order keys (the
    integer order `f64::total_cmp` uses) and bucketed by the high bits of `key − min`.
    The exponent sits in the high bits, so buckets adapt to skewed data. Value-range
    buckets were tried first and were no faster than sorting on log-normal data.
  - **Narrowing:** only buckets that contain a wanted rank are gathered, branch-free, and
    refined the same way, or finished by partial selection. Two adjacent ranks (a
    quantile's prev/next pair) cost one `select_nth_unstable` plus a minimum.
- **Identical because:**
  - The value at a given sorted rank does not depend on how it is found.
  - Equal keys are identical bit patterns.
  - `(x − mean) / scale` with a finite mean and a positive, finite scale preserves order,
    including −0.0 before +0.0. So the selected z-values are bit-for-bit those of the full
    sort.
- **Effect:** at T = 10,000, about 57–70 µs per column instead of about 265 µs for the
  sort (4–4.7×, `examples/select_bench.rs`).
- **Weak spot:** a dense bulk with about 1 % far outliers takes 121 µs (2.2×). The
  outliers stretch the first level's key range, so the bulk gets few buckets. Bounding
  the first level by sample quantiles would fix this, if real data shows the pattern.

### T3. Delta kernel: iterate rows in the outer loop

- **Where:** [src/delta.rs:190](src/delta.rs:190) (combine loop) and the per-column
  gather in `Reducer::reduce` ([src/delta.rs:259](src/delta.rs:259),
  [src/delta.rs:270](src/delta.rs:270)).
- **Now:** the loops go column by column over row-major T × b arrays, and index the
  (C, T) count matrix by residue code. Python sizes blocks to fit the cache, but the
  width is clamped at 16, so beyond about 15,600 taxa a block no longer fits: at
  200,000 taxa each block is 25.6 MB and is swept 16 times for the combine step and 16
  more for the reduce.
- **Change:** loop over taxa outside and columns inside; keep a (T, C) copy of the
  counts for the per-row lookup; write the result column-major so the reduce reads each
  column contiguously. Leave the BLAS call and the block boundaries exactly as they are.
- **Identical because:** each output element is computed by the same expression from the
  same inputs; only the order in which independent elements are visited changes. The
  degenerate-column mask is an OR, which is order-independent.
- **Effect:** memory traffic per block drops from about 1.2 GB to under 0.2 GB. This
  only shows above about 15,600 taxa, so the current benchmark does not see it.

### T4. Shuffle: remove the unpredictable branches, batch the counting

- **Where:** [src/rng.rs:127](src/rng.rs:127) `random_interval`,
  [src/rng.rs:150](src/rng.rs:150) `shuffle`,
  [src/chi2.rs:79](src/chi2.rs:79) `add_column_counts`.
- **Now:** each bounded draw runs a rejection loop whose branch is taken 0–50 % of the
  time, plus a second branch for the buffered 32-bit half. About 10 ns per cell, which
  points to branch mispredictions rather than arithmetic. The counts are then scattered
  into a 32 MB (C, T) array of `i64` once per column.
- **Change:**
  - Generate raw 32-bit draws into a small buffer and filter them without branches
    (store the candidate, advance the write index by the accept flag). The mask only
    changes when the bound crosses a power of two, so run one tight loop per power-of-two
    range.
  - Shuffle several columns, then add their counts in one sweep over the taxa; use
    `u32` counts and widen at the end.
- **Identical because:** the generator emits the same stream and the same draws are
  accepted; counts are integers.
- **Effect:** perhaps 2–3× on the null matrix and on the shuffle part of the
  delta-null. It also makes M1's second pass cheap.

### T5. Encoding: tiled, parallel transpose

- **Where:** [src/chi2.rs:27](src/chi2.rs:27) `from_rows`,
  [src/alignment.rs:152](src/alignment.rs:152) `detect_sequence_type`.
- **Now:** one thread writes the column-major matrix row by row, so consecutive writes
  are 200 KB apart: 8 × 10⁹ cache and TLB misses. The sequence-type scan is a second
  single-threaded pass.
- **Change:** transpose in cache-sized tiles, in parallel over column ranges; collect the
  set of residues present in the same pass that reads the file.
- **Identical because:** same codes in the same places.
- **Effect:** minutes → seconds at the target size.

### T6. Reuse kernel scratch buffers

- **Where:** `per_row_block` and `reduce` allocate and zero their work arrays on every
  block ([src/delta.rs:162](src/delta.rs:162)): two 25.6 MB arrays per block at 200,000
  taxa, across all threads.
- **Change:** per-thread scratch (rayon `map_init`), no zeroing where every element is
  overwritten.
- **Identical because:** buffers are fully written before they are read.
- **Effect:** probably 5–10 % of kernel time, and less page-fault contention at high
  thread counts.

### T7. Parallelism gaps

- **Few blocks at moderate taxon counts.** Prune iterations are parallel over column
  blocks only. With 1,000 taxa the block width is 250, so a 1,000-column alignment has 4
  blocks: the iteration stage went from 3.55 s to 1.04 s at 32 workers. Parallelizing
  the per-column reduce inside a block fixes this. Not an issue at 200,000 taxa (2,500
  blocks).
- **Delta-null with fewer permutations than threads.** Each permutation runs its kernel
  on one thread ([src/delta.rs:404](src/delta.rs:404)). Letting the blocks inside a
  permutation run in parallel uses the idle threads.
- **Serial summaries.** `NullState::new` sorts the pooled null twice
  ([src/permutation.rs:144](src/permutation.rs:144), `kth_smallest` at
  [src/npy.rs:102](src/npy.rs:102)) and takes two medians. One shared sorted copy serves
  the quantiles, the threshold and T1. Seconds, not minutes.
- **Identical because:** results are collected in index order in every case.

### T8. Build settings

- **Where:** [Cargo.toml](Cargo.toml) `[profile.release]`.
- **Change:** `lto = "fat"`, `codegen-units = 1`, and a newer CPU baseline
  (`-C target-cpu=x86-64-v3`, or runtime dispatch for the hot loops).
- **Identical because:** rustc does not reorder floating-point operations or fuse
  multiply-add on its own; wider instructions only change how independent elements are
  processed.
- **Effect:** a few percent to perhaps 20 % on the integer and element-wise loops. The
  v3 baseline needs AVX2 (any CPU that runs the Haswell/Zen OpenBLAS kernels has it).

## Looked at and left alone

- **Block width, chunking, BLAS call shapes.** These decide which OpenBLAS kernel and
  which numpy summation order apply, so they must stay as Python has them.
- **Incremental deltas across iterations.** Removing a column changes every taxon's
  counts, so every delta legitimately changes; no exact shortcut.
- **Skipping columns that cannot be the maximum.** Needs a bound that is sound in
  floating point; none is cheap.
- **Deduplicating identical columns.** A column's BLAS result can depend on its position
  inside the block, so duplicates in different positions may differ in the last bit.
- **Packing codes to 5 bits.** Saves 3 GB of 8 and slows every hot loop.

## Suggested order

1. T1, M2, M1 — small, contained, and they remove the hard blockers.
2. T2, T3, T6 — the kernel; verify against `tests/intermediates` and the goldens.
3. M3, M4, M5, T5 — the I/O path, done together because they share the storage change.
4. T4, T7, T8.

The existing parity suite (33 + 9 goldens, 15 intermediate dumps, the large tier and the
full-MSA reference) is the safety net for each step. Add a golden with more than 15,625
taxa so the large-block paths are exercised.
