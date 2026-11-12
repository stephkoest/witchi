# witchi-rust

A Rust port of [WitChi](https://github.com/stephkoest/witchi) v0.2.0, which tests and prunes
compositionally biased columns in multiple sequence alignments. On the same machine and with the
same OpenBLAS build, its output files are **byte-identical** to the Python reference. See
[PLAN.md](PLAN.md) for the design, the numerical hazards and how each one is handled.

## Build

```bash
cargo build --release
```

By default the binary links the OpenBLAS that numpy uses in the `witchi-0.2.0` conda environment.
That's required for bit-identical results. Point elsewhere with `OPENBLAS_LIB_DIR` (and
`OPENBLAS_LIB_NAME`, default `cblas`). For a portable build without OpenBLAS:

```bash
cargo build --release --no-default-features
```

The pure-Rust backend agrees with Python except, at most, in the last bits of a few Δ-scores.

## Usage

The command line matches the Python tool:

```bash
target/release/witchi prune --file alignment.fasta --pruning_algorithm wasserstein --top_n 10 --permutations 100
```

```bash
target/release/witchi test --file alignment.fasta --permutations 100 --create_output
```

Supported input formats are `fasta` and `nexus`.

The Rust port adds `--threads N` (default: the larger of the `--num_workers_*` values). Thread
count never changes results. `--num_workers_chisq` keeps its Python meaning, because Python's
column chunking depends on it. `WITCHI_TIMING=1` prints stage timings to stderr.

## Tests

```bash
cargo test --release
```

This runs the parity tests against Python fixtures: formatting, readers, RNG, intermediate
arrays, and the 33 `prune` and 9 `test` golden cases.

```bash
cargo test --release --test golden -- --ignored golden_large
```

That runs the 1000 × 6000 tier, after creating it with `scripts/make_golden.sh --large`.

Regenerate references from Python (`witchi-0.2.0` conda env):

```bash
scripts/make_golden.sh
```

```bash
scripts/make_golden.sh --test
```

```bash
scripts/regen_fixtures.sh
```
