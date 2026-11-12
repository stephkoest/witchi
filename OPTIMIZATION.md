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

Implemented: **T1, M2, M1** (PR #10), **T2** (PR #11), **T3 and T6** (PR #12), **T4** (PR #13),
**M3, M4, M5 and T5** (PR #14), then the main part of **T7** (PR #16). **T8** was measured
and left out (no gain). With them, all parity suites pass:
- the 33 + 5 + 9 goldens. The large tier now includes two cases at 20,000 taxa, above the
  15,625 where the block width is clamped at 16; both match Python byte for byte.
- the 15 intermediate cases
- the I/O, formatting and RNG fixtures
- the full-MSA run, byte-identical to Python

**Before vs after, and a profile near the target size:** see `bench/COMPARISON.md`.
- **Output:** on 18 cases from 100 to 200,000 taxa, the optimized binary gives the same
  output as the pre-optimization port (`84413af`).
- **Speed:** wall time is 2–5× lower at 10,000–50,000 taxa and 47× lower at 200,000 × 2,000.
- **Memory:** peak memory no longer grows with the alignment size multiplied by the threads.
- **Projection:** a profile at 200,000 × 10,000 with the defaults on 48 threads puts the
  target size at about 27 minutes and 18 GB. The delta-null (46 %) and the prune iterations
  (38 %) dominate.

Unit tests check, bit for bit:
- the sorted-null count against the elementwise count
- the streamed delta-null against the materialized one
- rank selection against a full sort (adversarial value sets, forced bucketing at small sizes)
- the Wasserstein reducer against the old sort-based reducer
- the column-major kernel and reducers against the old row-major ones (all algorithms,
  C- and F-ordered blocks, widths 1–37, row counts across the 128-row tile)
- the batched 32-bit draw stream and two-phase shuffle against numpy's draw-by-draw
  shuffle (sizes 0–70,001 around powers of two, many shuffles on one stream), and the
  batched permuted counts against column-by-column counting

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

### T3 and T6

None of the parity inputs exceeds 15,625 taxa, so these were also timed on a synthetic
51,776 × 2,000 input: the full MSA stacked four times, 2,000 evenly spaced columns,
20 permutations, cap 10, 32 threads. The output files are byte-identical to the T2
binary for all three algorithms.

| Input | Algorithm | Wall (s) | Delta-null (s) | Prune iterations (s) |
|---|---|---|---|---|
| 51,776 × 2,000 | wasserstein | 12.7 → 11.1 | 6.3 → 5.5 | 2.2 → 1.5 |
| 51,776 × 2,000 | squared | 16.2 → 11.2 | 9.0 → 5.1 | 3.0 → 1.4 |
| 51,776 × 2,000 | quartic | 17.0 → 10.3 | 10.0 → 5.0 | 2.8 → 1.3 |
| 10,000 × 10,000 (32 threads) | wasserstein | 29.4 → 26.7 | 15.7 → 14.4 | 6.8 → 5.6 |
| 10,000 × 10,000 (96 threads) | wasserstein | 21.6 → 19.1 | 10.9 → 9.5 | 5.5 → 4.5 |
| Full MSA 12,944 × 13,540 | wasserstein | 12.5 → 12.1 | 7.4 → 7.2 | 0.73 → 0.59 |

### T4

Same synthetic input and settings as T3/T6 (32 threads); every output file is byte-identical
to the previous binary for all three algorithms:

| Input | Algorithm | Wall (s) | Null matrix (s) | Delta-null (s) |
|---|---|---|---|---|
| 51,776 × 2,000 | wasserstein | 11.1 → 9.4 | 1.71 → 1.00 | 5.46 → 4.83 |
| 51,776 × 2,000 | squared | 11.2 → 8.8 | 2.18 → 0.72 | 5.12 → 4.14 |
| 51,776 × 2,000 | quartic | 10.3 → 7.3 | 1.66 → 1.04 | 4.99 → 3.06 |
| 10,000 × 10,000 (32 threads) | wasserstein | 26.7 → 24.8 | 4.32 → 2.85 | 14.4 → 12.8 |
| 10,000 × 10,000 (96 threads) | wasserstein | 19.1 → 17.1 | 2.53 → 1.52 | 9.50 → 8.41 |
| Full MSA 12,944 × 13,540 | wasserstein | 12.1 → 10.4 | 2.37 → 1.64 | 7.22 → 5.96 |

Peak RSS rises by about 5 % (1.49 → 1.57 GB on the synthetic input). Each permutation
in flight now holds a batch of 32 shuffled columns and a (T, C + 1) `u32` count table.
At 200,000 taxa that is about 23 MB per thread, against 32 MB for the `i64` counts it
already held.

### M3, M4, M5 and T5

The alignment now exists in memory once, as the encoded column-major matrix, from reading
to writing. Measured as before; all output files byte-identical (full MSA vs Python, the
others vs the previous binary and the original benchmark checksums):

| Input | Threads | Wall (s) | Peak RSS (MB) | Read (s) | Write outputs (s) |
|---|---:|---|---|---|---|
| Full MSA 12,944 × 13,540 | 32 | 10.4 → 10.4 | 1,055 → 528 | 0.52 → 0.28 | 0.70 → 0.23 |
| 10,000 × 10,000 | 32 | 24.8 → 23.4 | 855 → 467 | 0.38 → 0.20 | 1.29 → 0.45 |
| 10,000 × 10,000 | 96 | 17.1 → 16.3 | 1,533 → 1,144 | 0.37 → 0.24 | 1.32 → 0.43 |
| 10,000 × 100 | 32 | 1.78 → 0.94 | 419 → 374 | 0.01 → 0.01 | 0.94 → 0.20 |
| 51,776 × 2,000 (wasserstein) | 32 | 9.4 → 8.5 | 1,572 → 1,152 | | |
| 51,776 × 2,000 (squared) | 32 | 8.8 → 6.5 | 1,537 → 1,148 | | |

What is left of peak memory is the alignment (T × L bytes) plus per-thread working
memory, which grows with T. At the target size that is 8 GB plus about 50–60 MB per
thread at 200,000 taxa. Each permutation in flight holds:
- C × T counts and their f64 copies
- the kernel's two T × b blocks
- the counting batch and table

### T7 (segmented permutations) and T8

**T7.** Measured on the full MSA (20 permutations, 32 threads), alternating the binaries,
three runs each:
- **Delta-null:** 5.7–6.4 s → 4.5–4.9 s.
- **Null matrix:** 1.75 s → 1.3 s in two of three runs.

With 8 permutations on 8 threads (one segment), both stages are unchanged.

The 18 benchmark configurations on 10,000 taxa keep their checksums. Their timings were
taken while the machine's load average rose to 50, so they are only indicative: on
10,000 × 10,000 at 96 threads, delta-null 8.5 → 6.8 s and wall 16.3 → 14.5 s.

**T8.** Four builds were compared, at 4 threads (best of 3) and on the full MSA:
- the current settings
- `lto = "fat"` with `codegen-units = 1`
- that plus `-C target-cpu=x86-64-v3`
- that plus `-C target-cpu=native`

All gave the same output checksums. No variant was measurably faster: wasserstein took
3.02–3.11 s with every build, and on squared the current build was fastest. The hot paths
are OpenBLAS, memory-bound sweeps and code already written for the CPU's strengths. The
release profile is left as it is.

Remaining T7 items, not done: few blocks at moderate taxon counts, and serial summaries.
Neither matters at the target size (see T7).

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
- **Done** (`Encoded::from_alignment`, `chi2::Decoder`):
  - The matrix is encoded in place, with no copy.
  - Non-char-set bytes get distinct negative codes, -1 to at most -124.
  - The two histogram indexers clamp negative codes to -1. Every other kernel already
    tested `code >= 0`.

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
- **Done** (`fasta::read`, `AlignmentBuilder`), with the two-pass read:
  - **Pass 1** validates UTF-8 across the whole file, carrying incomplete characters
    across reads, and counts the records. A whole-file check keeps Python's error
    precedence.
  - **Pass 2** parses line by line with universal newlines. A lone `\r` ends a line,
    and `\r\n` counts as one terminator. Residues go straight into the column-major
    matrix (see T5).
  - **NEXUS** keeps its text parser and feeds the same builder.
  - The text parser (`fasta::parse`) is kept as the reference that a unit test compares
    the streaming reader against.

### M5. Output: stream the FASTA and the JSON

- **Where:** [src/output.rs:62](src/output.rs:62) (all pruned rows, then the whole FASTA
  in one buffer) and [src/output.rs:127](src/output.rs:127) (a JSON tree with one node
  per float: 20 million nodes for the permuted-null matrix).
- **Change:** write the FASTA through a buffered writer, transposing a tile of rows at a
  time from the column-major codes. Write the JSON incrementally, formatting each
  permutation's row of floats in parallel and emitting the rows in order.
- **Identical because:** same bytes in the same order.
- **Done:**
  - **FASTA:** groups of 64 taxa (one cache line of every column) are gathered from the
    codes through the `Decoder`, rendered in parallel and written in order.
  - **JSON:** written member by member. The (P, T) permuted null and the per-taxon lists
    go through `json_float_list`, which a unit test checks against the `Json` tree. The
    rows of the null are formatted in parallel.
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
- **Done:**
  - **Output:** the per-row block is written column-major (`per_row[k * T + t]`).
  - **Main loop:** it runs in tiles of 128 rows, so the strided reads of the BLAS
    product stay in cache while each column is written contiguously.
  - **Reducers:**
    - Wasserstein selects directly on a contiguous column slice, without copying it.
    - The pairwise sums read the slice.
    - The sequential sums run 8 columns side by side, each still adding rows in order
      t = 0, 1, ….
  - **Per-column counts:** an integer histogram, giving the same small-integer f64
    values as before.
  - **Unchanged:** the BLAS call and the block boundaries.
- **Identical because:**
  - Each element is computed by the same expression from the same inputs.
  - Each column's sum adds the same values in the same order as numpy's axis-0
    reduction of the row-major block.
  - The degenerate-column mask is an OR.
- **Effect:** at 51,776 taxa, about 2× on the kernel stages for squared and quartic, and
  about 1.4× on Wasserstein's pruning iterations. At 10,000 taxa, 8–12 % on wall time.
  See Status.

### T4. Shuffle: remove the unpredictable branches, batch the counting

- **Where:** [src/rng.rs:127](src/rng.rs:127) `random_interval`,
  [src/rng.rs:150](src/rng.rs:150) `shuffle`,
  [src/chi2.rs:79](src/chi2.rs:79) `add_column_counts`.
- **Now:** each bounded draw runs a rejection loop whose branch is taken 0–50 % of the
  time, plus a second branch for the buffered 32-bit half. About 10 ns per cell, which
  points to branch mispredictions rather than arithmetic. The counts are then scattered
  into a 32 MB (C, T) array of `i64` once per column.
- **Done** (`rng::U32Stream`, `permutation::permuted_counts`):
  - **Draw stream:** numpy's 32-bit draws (the low, then the high half of each 64-bit
    output) are generated 512 at a time.
  - **Two-phase shuffle:**
    - Phase 1 runs `random_interval`'s accept/reject process for every position without
      unpredictable branches. Within a power-of-two range the mask is constant, and a
      rejected draw's target is overwritten by the next one.
    - Phase 2 applies the swaps in numpy's order.
    - A first, fully branch-free version that swapped as it went was slower: 9.8 ns per
      cell against 7.4. Each rejected draw then waited on the swap's memory traffic.
  - **Batched counting:** 32 shuffled columns are added to a taxa-major (T, C + 1)
    `u32` table in one sequential sweep. Slot 0 takes codes outside the char set. The
    table is widened to (C, T) `i64` at the end. Both the null matrix and the
    delta-null's first pass use this.
- **Identical because:**
  - The stream is the same sequence of draws.
  - Phase 1 yields exactly the values `random_interval` returns, consuming the same
    draws.
  - The swaps are applied in the same order.
  - Counts are integers.
- **Effect:**
  - Shuffle: 7.4 → 4.1 ns per cell at T = 10,000, and 7.6 → 4.3 at T = 200,000
    (`examples/shuffle_bench.rs`).
  - Counting at T = 200,000: the per-column scatter cost 7.4–8.4 ns per cell, because
    the 32 MB `i64` table falls out of cache. The batched sweep removes that.
  - See Status for whole-run numbers.

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
- **Done:**
  - `AlignmentBuilder` buffers 256 rows and transposes them into the matrix in parallel,
    over blocks of 64 columns. Each block is a contiguous part of the matrix.
  - The set of bytes present is computed at the end, in parallel.
  - Encoding is an in-place, parallel table lookup.
  - Full MSA read: 0.28 s at 32 threads, 0.65 s on one thread.

### T6. Reuse kernel scratch buffers

- **Where:** `per_row_block` and `reduce` allocate and zero their work arrays on every
  block ([src/delta.rs:162](src/delta.rs:162)): two 25.6 MB arrays per block at 200,000
  taxa, across all threads.
- **Done:** a thread-local `KernelScratch` holds the block and reducer buffers, and the
  selection scratch, for the life of the thread. The BLAS output is still zero-filled
  before the call, as the freshly allocated array was.
- **Identical because:** every buffer is fully rewritten before it is read.
- **Effect:** measured together with T3 (see Status). Peak RSS is unchanged: each thread
  keeps about 2 × T × b × 8 bytes, which it already needed while processing a block.

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
- **Done: the permutation gap**, for the null matrix and the delta-null.
  - **Why it mattered:** each permutation ran on one thread. With P permutations on N
    threads, the last of ceil(P / N) rounds leaves threads idle. 100 permutations on
    96 threads ran in two rounds, the second with 4 permutations (about 52 % efficiency),
    and 20 permutations on 32 threads used only 20.
  - **How:** each permutation can now be split into column segments that run in
    parallel. A segment needs the random stream at its first column. A column's shuffle
    consumes a number of draws that depends only on the draws and the column length, so
    a draws-only pass finds each segment's position: `U32Stream::skip_shuffle`, about
    2.5 ns per cell against 4.4 ns for the shuffle. The segment then jumps its own PCG64
    there: `Pcg64::advance`, the standard O(log n) LCG jump-ahead (numpy's `advance`),
    and `U32Stream::at`.
  - **Combining:** counts are summed across segments. The delta-null's per-column
    deltas are concatenated in order, and Python's `max()` is taken over the whole list.
  - **Choosing the number of segments:** `segment_count` estimates the makespan
    `ceil(P * n / threads) * cost / n + skip` for n = 1..64 and picks the minimum, with
    per-cell costs of about 50 ns (delta-null), 6 ns (null matrix) and 2.5 ns (draws-only
    pass). With enough permutations per thread it stays at 1, the previous code path.
  - **Identical because:** every column is shuffled with exactly the draws of the
    sequential pass, and the results are combined in column order.
  - **Pitfall found:** the first version put the counting loop in a closure. It ran
    about 40 % slower, because sizes captured by reference were reloaded in the hot loop
    and blocked optimization. A plain function restores the old speed.
- **Not done:**
  - Parallelizing inside a block at moderate taxon counts: irrelevant at the target
    size, where there are 2,500 blocks.
  - A shared sorted null for the summaries: about seconds at P × T = 20 million.

### T8. Build settings

- **Where:** [Cargo.toml](Cargo.toml) `[profile.release]`.
- **Change:** `lto = "fat"`, `codegen-units = 1`, and a newer CPU baseline
  (`-C target-cpu=x86-64-v3`, or runtime dispatch for the hot loops).
- **Identical because:** rustc does not reorder floating-point operations or fuse
  multiply-add on its own; wider instructions only change how independent elements are
  processed.
- **Effect:** a few percent to perhaps 20 % on the integer and element-wise loops. The
  v3 baseline needs AVX2 (any CPU that runs the Haswell/Zen OpenBLAS kernels has it).
- **Measured: no gain**, and the output is identical (see Status). Not adopted.

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

1. T1, M2, M1 — small, contained, and they remove the hard blockers (done).
2. T2, T3, T6 — the kernel; verify against `tests/intermediates` and the goldens (done).
3. T4 (done).
4. M3, M4, M5, T5 — the I/O path, done together because they share the storage change (done).
5. T7 (permutation segmentation done), T8 (measured, no gain).

The existing parity suite (33 + 9 goldens, 15 intermediate dumps, the large tier and the
full-MSA reference) is the safety net for each step. The large tier now has two cases at
20,000 taxa (`gtdb_n20000_c2000`), so the clamped block-width paths are checked against
Python.
