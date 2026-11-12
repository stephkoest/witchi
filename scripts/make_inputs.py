#!/usr/bin/env python3
"""Generate the deterministic test alignments used by the golden tests.

Usage:
    python -I scripts/make_inputs.py [--large] SOURCE_MSA OUT_DIR

SOURCE_MSA is test_data/gtdb_r237_ar53_concatenated.masked.faa (single-line
FASTA, 12,944 taxa x 13,540 columns). Everything here is deterministic, so
re-running produces byte-identical files.

Small inputs (committed, tests/data/):
  gtdb_n20_c100      first 20 taxa, columns 0..100
  gtdb_n50_c500      first 50 taxa, columns 0..500
  gtdb_n100_c1000    every 129th taxon (diverse), columns 2000..3000
  gtdb_n50_c257      50 taxa, 257 columns: with T=50 the block width B is 256,
                     so the last column block has width 1 (numpy reductions
                     switch from sequential to pairwise for b == 1)
  gtdb_n150_c600     150 taxa x 600 columns: T*B*8 = 150*256*8 >= 256 KiB, so after
                     multi-column deletions numpy's temporary elision makes the
                     per-row block F-ordered (pairwise axis-0 sums; PLAN.md H16)
  gtdb_quirks        30 taxa x 400 columns with lowercase residues, multi-word
                     headers, wrapped/ragged sequence lines, blank lines and
                     non-standard residues (B, Z, *) -> "unexpected characters"
  dna_n30_c300       synthetic DNA with GC-biased taxa, gaps, N, an all-gap
                     column, an all-N column and duplicated columns (ties)
  dna_degenerate     DNA where 'T' occurs only in column 0 (once per taxon):
                     removing column 0 empties a character total and triggers
                     the direct per-column fallback path
Large inputs (--large, not committed, tests/golden_large/data/):
  gtdb_n1000_c6000   T*L >= 5e6, so --num_workers_chisq > 1 changes the column
                     chunking (PLAN.md hazard H15)
  gtdb_n20000_c2000  20,000 taxa x 2,000 evenly spaced columns. Above 15,625 taxa
                     the column-block width is clamped at 16 (a block no longer
                     fits Python's 2 MB budget) and rank selection uses its widest
                     bucketing. The source has 12,944 taxa: all are kept, and
                     7,056 synthetic taxa are added, each a copy of a source taxon
                     with 15 % of its columns replaced by the residue another
                     source taxon has in that column (so column compositions stay
                     realistic).
"""

import random
import sys
from pathlib import Path


def read_single_line_fasta(path, wanted_rows):
    """Return [(header, seq)] for the 0-based record indices in wanted_rows."""
    wanted = set(wanted_rows)
    out = {}
    with open(path) as fh:
        idx = -1
        header = None
        for line in fh:
            if line.startswith(">"):
                idx += 1
                header = line[1:].rstrip("\n")
                if idx > max(wanted):
                    break
            elif idx in wanted:
                out[idx] = (header, line.rstrip("\n"))
    return [out[i] for i in wanted_rows]


def write_fasta(path, records, wrap=None):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            if wrap is None:
                fh.write(seq + "\n")
            else:
                for i in range(0, len(seq), wrap):
                    fh.write(seq[i : i + wrap] + "\n")


def subset(src, rows, c0, c1):
    return [(h, s[c0:c1]) for h, s in read_single_line_fasta(src, rows)]


def make_quirks(src, out):
    rng = random.Random(7)
    recs = subset(src, list(range(200, 230)), 5000, 5400)
    with open(out, "w") as fh:
        for k, (h, s) in enumerate(recs):
            s = list(s)
            if k % 3 == 0:
                s = [c.lower() for c in s]
            if k % 7 == 1:  # sprinkle non-standard residues
                for _ in range(3):
                    s[rng.randrange(len(s))] = rng.choice("BZ*")
            s = "".join(s)
            title = h if k % 2 else f"{h} some description words\tk={k}"
            fh.write(f">{title}\n")
            if k % 5 == 4:
                fh.write("\n")  # blank line between header and sequence
            # ragged line lengths
            i = 0
            width = 50 + (k % 4) * 10
            while i < len(s):
                fh.write(s[i : i + width] + "\n")
                i += width
            if k % 6 == 5:
                fh.write("\n")


def make_dna(out):
    rng = random.Random(42)
    T, L = 30, 300
    biased = set(range(0, 30, 4))
    cols = []
    for j in range(L):
        col = []
        for t in range(T):
            if t in biased and j % 3 != 0:
                w = (0.1, 0.4, 0.4, 0.1)  # GC rich
            else:
                w = (0.3, 0.2, 0.2, 0.3)
            c = rng.choices("ACGT", weights=w)[0]
            r = rng.random()
            if r < 0.04:
                c = "-"
            elif r < 0.05:
                c = "N"
            col.append(c)
        cols.append(col)
    cols[17] = ["-"] * T  # all-gap column
    cols[90] = ["N"] * T  # all-N column
    for dst, src in ((200, 10), (201, 10), (250, 55)):  # exact duplicates -> tied deltas
        cols[dst] = list(cols[src])
    recs = [(f"dna_taxon_{t:02d}", "".join(cols[j][t] for j in range(L))) for t in range(T)]
    write_fasta(out, recs, wrap=70)


def make_dna_degenerate(out):
    rng = random.Random(3)
    T, L = 12, 60
    recs = []
    for t in range(T):
        body = ["A", "C", "G"] + [rng.choice("ACG") for _ in range(L - 4)]
        rng.shuffle(body)
        if t < 4:  # a little compositional bias
            body = [c if c != "A" or rng.random() < 0.3 else "G" for c in body]
            body[0] = "A"
            body[1] = "C"
        recs.append((f"deg{t}", "T" + "".join(body)))
    write_fasta(out, recs)


def spaced(n, total):
    return [i * total // n for i in range(n)]


def make_n20000(src, out):
    headers, seqs = [], []
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
            else:
                seqs.append(line.rstrip("\n"))
    cols = spaced(2000, len(seqs[0]))
    seqs = ["".join(s[j] for j in cols) for s in seqs]
    n_src, n_cols = len(seqs), len(cols)
    rng = random.Random(20000)
    recs = list(zip(headers, seqs))
    for k in range(20000 - n_src):
        base = (k * 7919) % n_src
        s = list(seqs[base])
        for j in rng.sample(range(n_cols), n_cols * 15 // 100):
            s[j] = seqs[rng.randrange(n_src)][j]
        recs.append((f"{headers[base]}_syn{k}", "".join(s)))
    write_fasta(out, recs)


def main():
    args = sys.argv[1:]
    large = "--large" in args
    args = [a for a in args if a != "--large"]
    if len(args) != 2:
        sys.exit(__doc__)
    src, out_dir = args
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if large:
        rows = [i * 12 for i in range(1000)]
        write_fasta(out / "gtdb_n1000_c6000.faa", subset(src, rows, 3000, 9000))
        make_n20000(src, out / "gtdb_n20000_c2000.faa")
        return
    write_fasta(out / "gtdb_n20_c100.faa", subset(src, list(range(20)), 0, 100))
    write_fasta(out / "gtdb_n50_c500.faa", subset(src, list(range(50)), 0, 500))
    write_fasta(
        out / "gtdb_n100_c1000.faa",
        subset(src, [i * 129 for i in range(100)], 2000, 3000),
    )
    write_fasta(
        out / "gtdb_n150_c600.faa",
        subset(src, [i * 80 + 3 for i in range(150)], 8000, 8600),
    )
    write_fasta(out / "gtdb_n50_c257.faa", subset(src, list(range(100, 150)), 7000, 7257))
    make_quirks(src, out / "gtdb_quirks.faa")
    make_dna(out / "dna_n30_c300.fasta")
    make_dna_degenerate(out / "dna_degenerate.fasta")


if __name__ == "__main__":
    main()
