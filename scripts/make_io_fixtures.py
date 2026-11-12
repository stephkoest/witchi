#!/usr/bin/env python3
"""Record how Python WitChi reads each IO test input.

Usage: python -I scripts/make_io_fixtures.py INPUT_DIR... --out EXPECTED_DIR

For every input (format by extension: .nex -> nexus, else fasta) this runs the
exact WitChi pipeline (AlignmentReader.run, SequenceTypeDetector.detect,
update_sequences, write_alignment) and writes EXPECTED_DIR/<name>.txt:

    OK <dna|protein|unexpected:<chars>|detect-error>
    <FASTA exactly as witchi would write the (unpruned) alignment>
or
    ERROR <ExceptionType>
"""

import contextlib
import io
import sys
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.alignment_reader import AlignmentReader
from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector


def describe(path, fmt):
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            alignment, array = AlignmentReader(str(path), fmt).run()
    except Exception as e:  # noqa: BLE001 - we record any failure
        return f"ERROR {type(e).__name__}\n"
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            is_dna, _ = SequenceTypeDetector.detect(alignment)
        seq_set = set()
        for record in alignment:
            seq_set.update(str(record.seq.upper()))
        seq_set -= {"N", "n", "-", "?", "X", "x", "."}
        if is_dna:
            kind = "dna"
        elif seq_set.issubset(set("ACDEFGHIKLMNPQRSTVWY")):
            kind = "protein"
        else:
            kind = "unexpected:" + "".join(sorted(seq_set))
    except ValueError:
        kind = "detect-error"
    pruner = AlignmentPruner.__new__(AlignmentPruner)
    records = pruner.update_sequences(alignment, array)
    buf = io.StringIO()
    AlignIO.write(MultipleSeqAlignment(records), buf, "fasta")
    return f"OK {kind}\n" + buf.getvalue()


def main():
    args = sys.argv[1:]
    out = Path(args[args.index("--out") + 1])
    inputs = [Path(a) for a in args[: args.index("--out")]]
    out.mkdir(parents=True, exist_ok=True)
    for d in inputs:
        for path in sorted(d.iterdir()):
            if path.name == "RUST_ONLY_ERRORS":
                continue
            fmt = "nexus" if path.suffix == ".nex" else "fasta"
            text = describe(path, fmt)
            (out / f"{path.name}.txt").write_bytes(text.encode())
            print(f"{path.name:28s} {text.splitlines()[0]}")


if __name__ == "__main__":
    main()
