#!/usr/bin/env python3
"""Generate the scaling-benchmark alignments (bench/data/, not committed).

Usage:
    python -I scripts/make_bench_inputs.py SOURCE_MSA OUT_DIR

For every n_taxa in (100, 1000, 10000) and n_cols in (100, 1000, 10000) writes
gtdb_t{n_taxa}_c{n_cols}.faa: taxa and columns are both evenly spaced across the
source MSA (test_data/gtdb_r237_ar53_concatenated.masked.faa, 12,944 x 13,540),
so each subset samples the whole taxonomy and all marker genes. Deterministic.
"""

import sys
from pathlib import Path

SIZES = (100, 1000, 10000)


def read_fasta(path):
    headers, seqs = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
                seqs.append([])
            else:
                seqs[-1].append(line.rstrip("\n"))
    return headers, ["".join(s) for s in seqs]


def spaced(n, total):
    return [i * total // n for i in range(n)]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out_dir = sys.argv[1:]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    headers, seqs = read_fasta(src)
    n_rows, n_cols_total = len(seqs), len(seqs[0])
    for t in SIZES:
        rows = spaced(t, n_rows)
        for c in SIZES:
            cols = spaced(c, n_cols_total)
            with open(out / f"gtdb_t{t}_c{c}.faa", "w") as fh:
                for r in rows:
                    s = seqs[r]
                    fh.write(f">{headers[r]}\n{''.join(s[j] for j in cols)}\n")
            print(f"gtdb_t{t}_c{c}.faa")


if __name__ == "__main__":
    main()
