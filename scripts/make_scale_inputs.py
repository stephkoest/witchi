#!/usr/bin/env python3
"""Generate large synthetic alignments for scale benchmarks (not committed).

Usage:
    python -I scripts/make_scale_inputs.py SOURCE_MSA OUT_DIR N_TAXA N_COLS [N_TAXA N_COLS ...]

Writes gtdb_syn_t{N_TAXA}_c{N_COLS}.faa for each pair. Needs numpy (run it in the
witchi-0.2.0 conda env).

SOURCE_MSA is test_data/gtdb_r237_ar53_concatenated.masked.faa (12,944 taxa x
13,540 columns). N_COLS columns are taken evenly spaced across the source. The
first min(N_TAXA, 12,944) taxa are the source taxa, evenly spaced. Each further
taxon k is a copy of source taxon (k * 7919) mod 12,944 with 15 % of its columns
replaced by the residue a random source taxon has in that column, so column
compositions stay realistic. Deterministic for a given numpy version.
"""

import sys
from pathlib import Path

import numpy as np


def read_source(path):
    headers, seqs = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
            else:
                seqs.append(line.rstrip("\n").encode("ascii"))
    return headers, np.frombuffer(b"".join(seqs), dtype=np.uint8).reshape(len(seqs), -1)


def generate(headers, msa, n_taxa, n_cols, path):
    n_src, total_cols = msa.shape
    cols = np.arange(n_cols) * total_cols // n_cols
    src = np.ascontiguousarray(msa[:, cols])
    rng = np.random.default_rng([n_taxa, n_cols])
    n_real = min(n_taxa, n_src)
    real = np.arange(n_real) * n_src // n_real
    chunk = 8192
    with open(path, "wb") as fh:
        for i in real:
            fh.write(b">" + headers[i].encode() + b"\n" + src[i].tobytes() + b"\n")
        for k0 in range(0, n_taxa - n_real, chunk):
            ks = np.arange(k0, min(k0 + chunk, n_taxa - n_real))
            base = (ks * 7919) % n_src
            block = src[base].copy()
            mask = rng.random(block.shape) < 0.15
            donors = rng.integers(0, n_src, size=block.shape)
            colidx = np.broadcast_to(np.arange(n_cols), block.shape)
            block[mask] = src[donors[mask], colidx[mask]]
            for k, b, row in zip(ks, base, block):
                fh.write(f">{headers[b]}_syn{k}\n".encode() + row.tobytes() + b"\n")


def main():
    args = sys.argv[1:]
    if len(args) < 4 or len(args) % 2:
        sys.exit(__doc__)
    src, out_dir = args[:2]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    headers, msa = read_source(src)
    pairs = args[2:]
    for t, c in zip(pairs[::2], pairs[1::2]):
        name = f"gtdb_syn_t{int(t)}_c{int(c)}.faa"
        generate(headers, msa, int(t), int(c), out / name)
        print(name)


if __name__ == "__main__":
    main()
