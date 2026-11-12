#!/usr/bin/env python3
"""Dump WitChi's internal arrays (bit-exact, as hex) for differential testing.

Usage (witchi-0.2.0 env):
    python -I scripts/dump_intermediates.py OUT_FILE --file F [--format fasta]
           [--pruning_algorithm A] [--permutations P] [--no-delta-null]

Each output line is ``name n hex1 hex2 ...`` with IEEE-754 bit patterns, so
Rust tests can compare without any decimal round-trip. Stages dumped:

  null_pooled, null_sums                 permutation null (P x T, P)
  z_mean, z_scale, z_quantiles, biased_threshold
  real_per_row_chi2, real_sum            the unpruned alignment
  initial_score, deltas                  first-iteration scores for the algorithm
  null_max_deltas                        delta-null distribution (if enabled)
"""

import argparse
import contextlib
import io
import struct

import numpy as np

from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.alignment_reader import AlignmentReader
from witchi.alignment_pruner.chi_square_calculator import (
    ChiSquareCalculator,
    _count_rows_from_int,
    _encode_alignment_int,
)
from witchi.alignment_pruner.permutation_test import PermutationTest
from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector


def hexes(values):
    arr = np.atleast_1d(np.asarray(values, dtype=np.float64)).ravel()
    return " ".join(f"{struct.unpack('<Q', struct.pack('<d', float(v)))[0]:016x}" for v in arr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--file", required=True)
    ap.add_argument("--format", default="fasta")
    ap.add_argument("--pruning_algorithm", default="wasserstein")
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--num_workers_chisq", type=int, default=1)
    ap.add_argument("--delta-null", action=argparse.BooleanOptionalAction, default=True)
    a = ap.parse_args()

    p = AlignmentPruner(
        file=a.file, format=a.format, permutations=a.permutations,
        pruning_algorithm=a.pruning_algorithm, delta_null=a.delta_null,
        num_workers_chisq=a.num_workers_chisq,
    )
    lines = []

    def dump(name, values):
        arr = np.atleast_1d(np.asarray(values, dtype=np.float64)).ravel()
        lines.append(f"{name} {arr.size} {hexes(arr)}".rstrip())

    with contextlib.redirect_stdout(io.StringIO()):
        alignment, arr = AlignmentReader(a.file, a.format).run()
        _, char_set = SequenceTypeDetector.detect(alignment)
        p.chi_square_calculator = ChiSquareCalculator(char_set, a.num_workers_chisq)
        p.permutation_test = PermutationTest(1, a.permutations)
        sums, pooled = p._compute_null_state(arr)
    dump("null_pooled", pooled)
    dump("null_sums", sums)
    dump("z_mean", p._null_z_mean)
    dump("z_scale", p._null_z_scale)
    dump("z_quantiles", p._null_z_quantiles)
    dump("biased_threshold", p._biased_taxa_threshold)

    p._alignment_int = _encode_alignment_int(arr, char_set)
    p._count_rows_raw = _count_rows_from_int(p._alignment_int, len(char_set))
    stats = p._calculate_per_row_stats(arr)
    dump("real_per_row_chi2", stats["per_row_chi2"])
    dump("real_sum", stats["sum"])
    initial, deltas = p.prune(arr, stats["expected_observed"], stats["count_rows"], pooled)
    dump("initial_score", initial)
    dump("deltas", [deltas[j] for j in range(len(deltas))])
    if p._null_max_deltas is not None:
        dump("null_max_deltas", p._null_max_deltas)

    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
