# Developing WitChi

WitChi 0.3 is a Rust port of the Python WitChi v0.2.0. On the same machine and with the same
OpenBLAS build, its output files are **byte-identical** to the Python reference, apart from the
version string: the port reports its own version (for example 0.3.1) in the console banner and in
the JSON `witchi_version` field, where Python 0.2.0 writes `0.2.0`.

- [PLAN.md](PLAN.md): the design, the numerical hazards and how each one is handled.
- [OPTIMIZATION.md](OPTIMIZATION.md): the performance work, all of which keeps the output
  byte-identical.
- [bench/](bench/): the scaling benchmark and the before/after comparison.

## Build

```bash
conda env create -f environment.yml
conda activate witchi
cargo build --release
```

By default the binary links `libcblas.so` from the active conda environment
(`$CONDA_PREFIX/lib`). `environment.yml` pins the OpenBLAS build (0.3.34) that Python WitChi
0.2.0 used for the reference outputs, which is required for bit-identical results. Point
elsewhere with `OPENBLAS_LIB_DIR` (and `OPENBLAS_LIB_NAME`, default `cblas`); the binary embeds
an rpath to that directory. OpenBLAS picks its kernels by CPU, so bit-identity also assumes the
same CPU family as the reference machine (the references were made on AMD EPYC 7402).

For a portable build without OpenBLAS:

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

That runs the large tier (1,000 × 6,000 and 20,000 × 2,000), after creating its inputs with
`scripts/make_inputs.py --large` and its references with `scripts/make_golden.sh --large`.

Regenerate references from Python (a `witchi-0.2.0` conda env with Python WitChi 0.2.0 installed):

```bash
scripts/make_golden.sh
```

```bash
scripts/make_golden.sh --test
```

```bash
scripts/regen_fixtures.sh
```
