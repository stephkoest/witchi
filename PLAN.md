# WitChi 0.2.0 → Rust port plan (`prune` first)

Reference implementation: `../witchi` (git `6c99225`, v0.2.0), installed in conda env
`witchi-0.2.0` (Python 3.11.18, numpy 2.4.6, Biopython 1.88, OpenBLAS 0.3.34 pthreads).
Dev machine: AMD EPYC 7402 (Zen 2, AVX2 + FMA), 96 threads, rustc 1.91.

Goal: `witchi-rs prune ...` writes the same four output files as Python, **byte for byte**.

---

## 1. What `prune` actually does (call graph)

```
main → AlignmentPruner.run
  AlignmentReader.run          Bio.AlignIO.read(fasta) → (T×L) '<U1' array, np.char.upper
  SequenceTypeDetector.detect  DNA "ACGT" vs AA "ACDEFGHIKLMNPQRSTVWY"
  _compute_null_state
    PermutationTest.compute_null   P permutations, seed 12345+i, permute each column
                                   → per-taxon χ² matrix (P×T); sums, pooled null
    null_z mean / MAD scale, K=min(200, P*T) quantiles of null Z
    _compute_biased_taxa_threshold  order statistic of the pooled null
    _compute_null_deltas (delta-null on)  P more permutations, seed 12345+P+i,
                                   max per-column Δ-score with the active algorithm
  recursive_prune               main loop (see §4)
  outputs                       *_pruned.fasta, *_pruned.tsv, *_pruned_scores.tsv,
                                *_pruned_score_dict.json   (prefix = <stem>_<algo>_s<top_n>)
```

Δ-scores per column (`_reduce_chi2_per_col_removed`) use a closed-form "remove column j"
identity evaluated in column blocks of width `B = max(16, min(256, 2_000_000 // (8*T)))`,
with a BLAS matmul `U.T @ inv_gf` per block and a direct fallback for degenerate columns.
Reducers: sum (squared), sum of squares (quartic), Wasserstein-1 on Z quantiles (wasserstein).

## 2. Byte-identity hazards and how each is handled

Each hazard below was checked against the installed environment unless noted otherwise.

| # | Hazard | Approach | Status |
|---|--------|----------|--------|
| H1 | `np.random.default_rng(seed)` + `apply_along_axis(rng.permutation, 0, a)` | Implement SeedSequence (hashmix pool of 4 u32) → PCG64 (128-bit LCG, XSL-RR) with numpy's buffered `next_uint32` → `random_interval` (mask rejection, 32-bit draws when bound ≤ 2³²−1) → Fisher–Yates `for i in (1..n).rev()` per column, one RNG across all columns | **Verified**: a pure-Python version matches numpy exactly (`probe.py`) |
| H2 | `U.T @ inv_gf` goes to OpenBLAS `dgemm`. The result is an FMA chain for interior tiles, but **edge tiles and the vector cases differ** (gemv/dot) | Call the **same OpenBLAS** through FFI, using numpy's matmul dispatch: `dot` if T==1&&b==1, `gemv` if T==1 or b==1, otherwise `dgemm(TransA=Trans, ...)` with the same strides. Single-threaded (Python uses `threadpool_limits(1)`). | **Verified hazard**: 59–100% of elements match a plain FMA chain depending on shape, so a pure-Rust kernel is *not* exact |
| H3 | numpy summation order: contiguous-axis reductions use **pairwise sum** (8 accumulators, block 128). Non-contiguous-axis reductions are **sequential**. When a block has width b==1 the axis-0 reductions become contiguous, so they switch to pairwise. | `npy::pairwise_sum` and `npy::seq_sum`. Each reduction call site is labeled with the variant numpy uses for that array layout (table in `src/npy.rs` docs). | Pairwise **verified** bit-exact against `np.sum` up to n=123457 |
| H4 | `np.quantile(..., method="linear")`: virtual index `(n-1)*q`, floor/gamma, `_lerp` with the `t >= 0.5` branch (`b - d*(1-t)`) | Port numpy 2.4.6 `_quantile` / `_lerp` line by line. Use sort or `select_nth` to get order statistics (the values are identical). | To do: copy from the installed numpy source |
| H5 | `np.linspace(0,1,K+2)[1:-1]`: `arange*step + start` and the last element forced | Port the formula exactly | To do |
| H6 | `np.median` (partition, then `(a+b)/2` for even n); `np.mean` = pairwise_sum/n; `np.percentile` (stdout only) | Port | To do |
| H7 | Float text: `str(np.float64)` / Python `repr` (shortest round-trip, sci notation when exp < −4 or ≥ 16, `1e-05`, `1e+16`, `inf`, `nan`); `json.dump` writes `NaN` / `Infinity` | `pyfmt::repr_f64` built on Rust's shortest `{:e}` digits, plus a hand-rolled JSON writer (indent=4, `", "`/`": "` separators, **no trailing newline**) | To do. Tests will fuzz random bit patterns against Python |
| H8 | `:.2f/.3f/.6f` formatting (stop_reason goes into the JSON) | Rust `{:.N}` rounds half-to-even on the exact value, like Python. Special-case inf/nan | Verify in tests |
| H9 | `csv` module: **`\r\n` line endings**, QUOTE_MINIMAL (quote fields containing `\t`, `"`, `\r`, `\n`), `None`→`""` | Small csv writer | Observed in the outputs |
| H10 | Biopython FASTA write: `>title` (id/description rule), sequence wrapped at 60 | Port `as_fasta` | Observed |
| H11 | Biopython FASTA read: id = first word of title, description = full title, sequence lines concatenated with whitespace removed, unequal lengths → error | Port `SimpleFastaParser` semantics | To do |
| H12 | Sort order: `sorted(..., key=delta, reverse=True)` is **stable**, so ties keep the lower column index first; `max(deltas.values())` | Stable sort by descending delta; the same `>` scan for max | – |
| H13 | `prune_dict` is an insertion-ordered dict; rollback deletes the last n inserted keys | `indexmap::IndexMap` + `shift_remove` | – |
| H14 | `np.char.upper` and the `cp & 0xFF` codepoint LUT in `_encode_alignment_int` alias non-ASCII chars (e.g. U+0141 behaves like 'A') | Read as bytes and reject non-ASCII with an error | Decided (§6) |
| H15 | Python's chunking with `--num_workers_chisq` > 1 (when T·L ≥ 5e6) changes the **block boundaries**, which changes BLAS edge tiles and b==1 blocks. Output can therefore depend on the worker count. | Reproduce `np.linspace(0, L, n_chunks+1, dtype=int)` chunk bounds and the per-chunk block loop exactly. Parallelise *across* chunks/blocks without changing their boundaries. | Must test with a large subset |

Python itself is only reproducible on the same CPU family and the same OpenBLAS kernel
(H2). "Byte identical" therefore means: same machine class, same OpenBLAS build.

## 3. Crate layout

```
witchi-rust/
  Cargo.toml            bin "witchi" (name TBD), deps: clap, indexmap, rayon, anyhow
  build.rs              link libopenblas (conda env lib dir via env var / pkg-config)
  src/main.rs           CLI: prune (and later test); argparse-compatible flag names
  src/fasta.rs          Biopython-compatible reader/writer
  src/seqtype.rs        SequenceTypeDetector
  src/rng.rs            SeedSequence, PCG64, random_interval, permute_columns
  src/npy.rs            numpy-exact numerics: pairwise_sum, seq_sum, mean, median,
                        kth/partition, quantile(linear), linspace, percentile
  src/blas.rs           numpy-matmul dispatch → cblas_ddot/dgemv/dgemm (FFI);
                        optional pure-Rust FMA fallback behind a feature flag
  src/chi2.rs           ChiSquareCalculator + block kernel + reducers
  src/permutation.rs    PermutationTest null (+ calc_empirical_pvalue)
  src/pruner.rs         AlignmentPruner state machine (§4)
  src/output.rs         4 output writers, csv/json writers
  src/pyfmt.rs          Python repr / format-spec helpers
  tests/golden/         small inputs + Python reference outputs
  tests/golden.rs       runs the binary, then byte-compares against the references
  scripts/make_golden.sh  regenerates the references from the conda env
  scripts/dump_intermediates.py  numpy dumps (hex floats) for differential debugging
```

Data representation: alignment as `Vec<u8>` row-major T×L (ASCII). Int-encoded as `i8`
(−1 = not in char set). Counts as `i64` (C×T, row-major like numpy). Every float array
keeps numpy's memory layout, so the reduction rules (H3) can be read straight off the
Python code.

## 4. Pruner state machine (port of `recursive_prune`)

Port literally, keeping every quirk:
- loop guard `removed + len(top_n_indices) < max_residue_pruned`
- iteration 0 records `before_real` / `_before_stats` (computed against the **original** null)
- prune_dict rows are written in the *next* iteration (global χ² after removal, biased
  count against the *current* null threshold, `null_pvals_for_batch.get(col)`)
- the `+1` count fudge is re-evaluated every iteration on cached raw counts; raw counts are
  decremented incrementally
- stopping: alignment p > 0.95 (delta-null on) or > 0.05 (off)
- delta-null walk: `p = #(null_max >= delta)/P`, break at p > 0.05
- touchdown rollback (once): drop the last batch from prune_dict, rebuild from the original
  alignment, `top_n = max(1, initial//10)`
- one-shot null recompute (`null_recompute_budget=1`), reusing seeds 12345+i, which also
  overwrites the biased threshold, the Z params and the delta-null
- natural-exit (cap) bookkeeping block, stop_reason strings verbatim
- `self.top_n` reset to the initial value afterwards. The output filename uses `self.top_n`
  *after* the reset, i.e. the initial value.

## 5. Phases and acceptance tests

0. **Harness** – `cargo init`; `make_golden.sh` builds subsets from
   `test_data/gtdb_r237_ar53_concatenated.masked.faa` (12,944 taxa × 13,540 cols, AA) and runs
   Python over a scenario matrix. Proposed subsets:
   `n20_c100`, `n50_c500`, `n100_c1000`, plus a synthetic DNA one and edge cases (all-gap
   column, taxon with a single residue → degenerate-column fallback, L ≡ 1 mod B → b==1
   block, lowercase input, multi-word FASTA headers).
   Scenarios: {squared, quartic, wasserstein} × {top_n 1, 5, 20 (touchdown)} ×
   {delta-null, no-delta-null} × {cap hit, natural stop}, with permutations 10–100.
1. **IO + formatting**: fasta round-trip, detector, `pyfmt` fuzz tests against Python.
   1b. **NEXUS reader**: parity with `AlignIO.read(..., "nexus")` on `example.nex` and
   on generated interleaved/quoted-label variants.
2. **RNG**: permuted alignments bit-identical to numpy for several seeds/shapes.
3. **Null**: pooled null matrix, sums, threshold, Z params, quantiles bit-identical
   (compare hex-float dumps from `dump_intermediates.py`).
4. **Δ-scores**: per-column deltas for all three algorithms, including fallback columns and
   b==1 blocks, bit-identical. Includes `blas.rs`.
5. **Pruner loop + outputs**: all golden scenarios byte-identical (`cmp`).
6. **Console output**: informative and close to Python's (not byte-checked).
7. **Scale + performance**: rayon over permutations and over column chunks; run
   medium (1k × 5k) and then the full dataset against Python; check H15 with
   `--num_workers_chisq` > 1; measure speed-up and memory.
8. **`test` subcommand** (later): reuses rng/null/npy/pyfmt; only the
   `_scores.tsv` writer and the console report are new.

## 5a. Phase 0 status (done)

- The installed reference package is tag v0.2.0 (`22bd7da`). It is identical to
  `../witchi/src` (HEAD `6c99225` only adds docs).
- `scripts/make_inputs.py`: deterministic subsets and synthetic edge-case inputs → `tests/data/`.
- `tests/cases.tsv`: 29 scenarios. `scripts/make_golden.sh` runs Python and writes
  `tests/golden/<case>/` (4 output files + informational `stdout.txt`) and `PROVENANCE-*.txt`.
- `tests/golden.rs`: byte-compares the Rust outputs (`#[ignore]` until phase 5).
  `WITCHI_GOLDEN_BIN=<python witchi>` re-runs Python through the harness, and all 29 cases
  reproduce, so the references are deterministic.
- Path coverage (`scripts/path_coverage.py` + stdout markers):
  - alignment-p stop, delta-not-significant stop, cap stop
  - touchdown via both alignment overshoot and delta-null partial batch
  - non-rollback partial batch (`n50_wa_t20`), null recompute, `--no-delta-null`
  - degenerate-column fallback (`deg_*`), b == 1 blocks (`c257_*`, `dna_sq_t1`)
  - K < 200 (`nex_*`: P·T = 100), unexpected residues (`quirks_*`), NEXUS with `matchchar=.`
- Large tier (`tests/cases_large.tsv`, 1000 × 6000, refs not committed) covers H15 chunking.
  Python's output was byte-identical for `--num_workers_chisq` 1 vs 4 here, even though
  chunk bounds stop being multiples of B once columns are pruned. H15 may be harmless in
  practice; the port still reproduces the chunk/block boundaries because it costs nothing.
- **Python crash found**: `--max_residue_pruned 0` raises `UnboundLocalError` (`per_row_chi2`
  is unset in `recursive_prune`). Under §6 the Rust port exits with a clear error instead.

## 6. Decisions (agreed 2026-11-12)

- **BLAS**: both backends. The OpenBLAS FFI backend is the default and gives exact parity.
  A pure-Rust FMA backend sits behind a cargo feature (e.g. `--no-default-features
  --features pure-rust`) for portable builds. Golden tests are exact for the OpenBLAS
  backend and use a ulp tolerance for the pure-Rust one.
- **Byte-identity scope**: the four output files only. Console output only needs to be
  informative (keeping it close to Python's is still useful for side-by-side debugging).
  Phase 6 becomes "reasonable console output".
- **Input formats**: FASTA + NEXUS. NEXUS must reproduce what `Bio.AlignIO.read(f, "nexus")`
  produces: taxon ids, the sequence after Biopython's missing/gap handling, and record
  order. Reference input: `../witchi/tests/data/example.nex`. This adds a phase **1b**:
  port the needed subset of `Bio.Nexus` (DATA/CHARACTERS block, MATRIX, interleave,
  quoted labels, comments, `format gap=/missing=`). Anything outside that subset is a
  clear error.
- **Pathological input**: exact parity for valid ASCII alignments. Clear error messages
  for non-ASCII bytes, unequal lengths, empty or zero-column alignments, and unknown
  `--format` / `--pruning_algorithm` values. Python quirks for those cases are not
  reproduced.

## 7. Still open (non-blocking; defaults in brackets)

- Binary name [`witchi`, same CLI flags as Python] and `git init` of this repo [yes].
- Where OpenBLAS is found at build time [`OPENBLAS_LIB_DIR`, defaulting to the
  `witchi-0.2.0` conda env's `lib/`, with an rpath so the binary runs without activating it].
