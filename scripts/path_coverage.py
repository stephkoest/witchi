#!/usr/bin/env python3
"""Report which internal numeric paths of Python WitChi a prune run exercises.

Usage (inside the witchi-0.2.0 conda env, from a scratch directory):
    python -I path_coverage.py <witchi prune args...>

stdout does not show these paths, so they are counted by wrapping the internals:
  * degenerate columns that take the direct per-column fallback (zero_mask)
  * column blocks of width b == 1 (numpy axis-0 reductions turn pairwise)
  * BLAS call shapes (T, b) passed to U.T @ inv_gf
  * the number of chunks used by _reduce_chi2_per_col_removed (H15)
"""

import collections
import sys

import numpy as np

import witchi.alignment_pruner.chi_square_calculator as csc
from witchi.witchi import main

counts = collections.Counter()
orig_block = csc._per_row_chi2_block


def block(cr, alignment_int, U, tcb, nb, ctb, j_start, j_end):
    b = j_end - j_start
    counts[f"block_width_{'1' if b == 1 else 'gt1'}"] += 1
    out, zero_mask = orig_block(cr, alignment_int, U, tcb, nb, ctb, j_start, j_end)
    counts["degenerate_columns"] += int(np.count_nonzero(zero_mask))
    return out, zero_mask


csc._per_row_chi2_block = block

orig_parallel = csc.Parallel


class CountingParallel(orig_parallel):
    def __init__(self, n_jobs=None, **kw):
        counts[f"chi2_parallel_chunks_{n_jobs}"] += 1
        super().__init__(n_jobs=n_jobs, **kw)


csc.Parallel = CountingParallel

sys.argv = ["witchi", "prune"] + sys.argv[1:]
try:
    main()
finally:
    print("PATH-COVERAGE", dict(sorted(counts.items())), file=sys.stderr)
