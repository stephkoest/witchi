"""
MSA Isolation Stratification
=============================

Assigns aligned sequences to strata based on evolutionary isolation
(min nearest-neighbour distance under gap-aware Hamming), for use in
stratified permutation testing of compositional bias.

Pipeline
--------
1. **Pairwise distance matrix** — gap-aware Hamming distance between
   every pair of sequences.                           O(N^2 * L)
2. **Isolation metric** — distance to the single nearest neighbour
   per taxon (min off-diagonal row-wise).  Direct proxy for terminal
   branch length: under stationary models, dist(taxon, sister) is
   approximately the sum of the two terminal branches.   O(N^2)
3. **Unsupervised binning** — Fisher-Jenks natural breaks on the 1-D
   isolation metric, with Calinski-Harabasz index for automatic K
   selection.  No thresholds.                         O(N*K_max^2)

Total: **O(N^2 * L)** — quadratic in N, linear in alignment length.

For witchi's typical N range (a few hundred to a few thousand), this is
a few seconds at most.  Above N ~ 5000, an LSH-Hamming or HNSW index
could replace the brute pairwise pass without changing the metric.

Requirements: numpy only.
"""

from __future__ import annotations

from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Encoding helpers
# ---------------------------------------------------------------------------

_GAP_ORDS = np.array([ord(c) for c in "-.*X"], dtype=np.uint8)


def _encode_msa(msa: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Encode full MSA once -> (uint8 array (N,L), bool residue mask (N,L))."""
    n = len(msa)
    L = len(msa[0])
    enc = np.empty((n, L), dtype=np.uint8)
    for i, seq in enumerate(msa):
        enc[i] = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    is_res = ~np.isin(enc, _GAP_ORDS)
    return enc, is_res


# ---------------------------------------------------------------------------
# Phase 1: Pairwise gap-aware Hamming  — O(N^2 * L), memory O(N^2)
# ---------------------------------------------------------------------------


def _hamming_seed_to_all(
    enc: np.ndarray, is_res: np.ndarray, seed_idx: int
) -> np.ndarray:
    """Fractional Hamming from one row to all rows (fully vectorised).

    Parameters
    ----------
    enc : (N, L) uint8 — pre-encoded MSA
    is_res : (N, L) bool — residue mask (True = not a gap)
    seed_idx : index of the reference row
    """
    shared = is_res & is_res[seed_idx]  # (N, L) bool
    n_shared = shared.sum(axis=1).astype(np.float64)  # (N,)
    mismatches = (
        ((enc != enc[seed_idx]) & shared).sum(axis=1).astype(np.float64)  # (N,)
    )
    return np.where(n_shared > 0, mismatches / n_shared, 1.0)


def _pairwise_hamming(enc: np.ndarray, is_res: np.ndarray) -> np.ndarray:
    """Full N x N gap-aware Hamming distance matrix.

    Each cell D[i, j] = (mismatch positions where both i and j have a
    residue) / (positions where both have a residue).  Diagonal = 0.
    """
    n = enc.shape[0]
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        D[i] = _hamming_seed_to_all(enc, is_res, i)
    return D


# ---------------------------------------------------------------------------
# Phase 2: Isolation metric — min off-diagonal per row  — O(N^2)
# ---------------------------------------------------------------------------


def _min_nn_distances(D: np.ndarray) -> np.ndarray:
    """Per-row min off-diagonal distance — distance to the nearest neighbour.

    Parameters
    ----------
    D : (N, N) symmetric distance matrix with zero diagonal.
    """
    masked = D.copy()
    np.fill_diagonal(masked, np.inf)
    return masked.min(axis=1)


# ---------------------------------------------------------------------------
# Phase 3: Fisher-Jenks natural breaks + Calinski-Harabasz  — O(N*K^2)
# ---------------------------------------------------------------------------


def _fisher_jenks(values: np.ndarray, k: int) -> np.ndarray:
    """Fisher-Jenks optimal 1-D partition into k classes (vectorised DP).

    Parameters
    ----------
    values : (N,) sorted array of floats
    k : number of classes (>= 2)

    Returns
    -------
    labels : (N,) array of class indices 0..k-1
    """
    n = len(values)
    assert k >= 2 and k <= n

    # Precompute cumulative sums for O(1) SSD of any contiguous slice
    cum = np.empty(n + 1, dtype=np.float64)
    cum2 = np.empty(n + 1, dtype=np.float64)
    cum[0] = 0.0
    cum2[0] = 0.0
    np.cumsum(values, out=cum[1:])
    np.cumsum(values * values, out=cum2[1:])

    def _ssd_vec(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        """Vectorised SSD for arrays of (lo, hi) pairs."""
        cnt = (hi - lo + 1).astype(np.float64)
        s = cum[hi + 1] - cum[lo]
        s2 = cum2[hi + 1] - cum2[lo]
        return s2 - s * s / cnt

    # dp_prev[i] = min SSD for values[0..i] using 1 class
    arange_n = np.arange(n)
    dp_prev = _ssd_vec(np.zeros(n, dtype=np.intp), arange_n)

    back = np.zeros((k, n), dtype=np.intp)

    for j in range(1, k):
        dp_curr = np.full(n, np.inf)
        back_j = np.full(n, j, dtype=np.intp)

        # For each position i (from j to n-1), find optimal split m in [j..i]
        # cost(m, i) = dp_prev[m-1] + ssd(m, i)
        # Vectorise over m for each i: precompute ssd(m, i) for all valid m
        # Process column-wise: for each candidate split m, compute cost for all i >= m
        for m in range(j, n):
            # ssd(m, i) for all i from m to n-1
            i_range = arange_n[m:]
            m_arr = np.full(len(i_range), m, dtype=np.intp)
            costs = dp_prev[m - 1] + _ssd_vec(m_arr, i_range)
            mask = costs < dp_curr[m:]
            dp_curr[m:][mask] = costs[mask]
            back_j[m:][mask] = m
        dp_prev = dp_curr
        back[j] = back_j

    # Backtrack to recover class labels
    labels = np.empty(n, dtype=np.intp)
    end = n
    for j in range(k - 1, 0, -1):
        start = int(back[j, end - 1])
        labels[start:end] = j
        end = start
    labels[0:end] = 0

    return labels


def _calinski_harabasz(values: np.ndarray, labels: np.ndarray) -> float:
    """Calinski-Harabasz index for a 1-D partition.

    CH = [B / (K-1)] / [W / (N-K)]

    where B = between-cluster SS, W = within-cluster SS.
    Higher is better.
    """
    n = len(values)
    k = int(labels.max()) + 1
    if k < 2 or n <= k:
        return 0.0

    grand_mean = values.mean()
    B = 0.0
    W = 0.0
    for c in range(k):
        mask = labels == c
        nc = mask.sum()
        if nc == 0:
            continue
        cluster_mean = values[mask].mean()
        B += nc * (cluster_mean - grand_mean) ** 2
        W += ((values[mask] - cluster_mean) ** 2).sum()

    if W == 0.0:
        return float("inf")
    return (B / (k - 1)) / (W / (n - k))


def _auto_strata(
    isolation: np.ndarray,
    min_stratum_size: int = 2,
    max_clusters: Optional[int] = None,
) -> np.ndarray:
    """Assign taxa to strata by isolation metric using Jenks + CH.

    Parameters
    ----------
    isolation : (N,) mean k-NN distances
    min_stratum_size : minimum taxa per stratum (>= 2)
    max_clusters : upper bound on K (default: min(20, N // min_stratum_size))

    Returns
    -------
    bin_ids : (N,) array of stratum indices 0..K-1
    """
    n = len(isolation)
    if max_clusters is None:
        max_clusters = min(20, n // min_stratum_size)
    max_clusters = max(1, min(max_clusters, n // min_stratum_size))

    if max_clusters < 2:
        return np.zeros(n, dtype=np.intp)

    # Sort isolation values, keep track of original indices
    order = np.argsort(isolation)
    sorted_vals = isolation[order]

    best_ch = -1.0
    best_labels = np.zeros(n, dtype=np.intp)

    for k in range(2, max_clusters + 1):
        labels_sorted = _fisher_jenks(sorted_vals, k)

        # Check min_stratum_size constraint
        valid = True
        for c in range(k):
            if (labels_sorted == c).sum() < min_stratum_size:
                valid = False
                break
        if not valid:
            continue

        ch = _calinski_harabasz(sorted_vals, labels_sorted)
        if ch > best_ch:
            best_ch = ch
            best_labels = labels_sorted.copy()

    # Map back to original taxon order
    result = np.empty(n, dtype=np.intp)
    result[order] = best_labels
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════


def msa_strata(
    msa: list[str],
    names: Optional[list[str]] = None,
    min_stratum_size: int = 2,
    max_clusters: Optional[int] = None,
    n_neighbours: Optional[int] = None,  # deprecated; kept for API compat
    n_seeds: Optional[int] = None,  # deprecated; kept for API compat
    return_diagnostics: bool = False,
) -> dict[str, int] | tuple[dict[str, int], dict]:
    """Assign each sequence to a stratum for stratified permutation testing.

    Strata group taxa by evolutionary isolation (min nearest-neighbour
    distance under gap-aware Hamming), a direct proxy for terminal branch
    length.  Within-stratum permutations then preserve the per-taxon
    drift profile imposed by branch-length heterogeneity.  The number of
    strata is selected automatically via Calinski-Harabasz on Fisher-Jenks
    natural breaks.

    Parameters
    ----------
    msa : list of str
        Aligned sequences (equal length, gaps allowed).
    names : list of str, optional
        Sequence identifiers.
    min_stratum_size : int
        Minimum sequences per stratum (>= 2 for permutation; default 2).
    max_clusters : int, optional
        Upper bound on number of strata.  Default: min(20, N // min_stratum_size).
    n_neighbours : int, optional
        Deprecated.  Previously controlled k-NN graph connectivity over a
        seed projection; the projection is no longer used.  Accepted but
        ignored.
    n_seeds : int, optional
        Deprecated.  Previously controlled seed-projection dimensionality;
        the projection is no longer used.  Accepted but ignored.
    return_diagnostics : bool
        If True, return (strata_dict, diagnostics_dict) instead of just
        strata_dict.  Diagnostics include the isolation metric and the
        number of natural strata (min_stratum_size=2) for comparison.

    Returns
    -------
    dict mapping sequence name (str) -> stratum index (int)
        — or (dict, diagnostics) if return_diagnostics=True
    """
    del n_neighbours, n_seeds  # deprecated; isolation is now exact min-NN

    n = len(msa)
    if names is None:
        names = [str(i) for i in range(n)]
    assert len(names) == n
    assert len({len(s) for s in msa}) == 1, "Sequences must be aligned (same length)"
    assert min_stratum_size >= 2, "min_stratum_size must be >= 2"

    if n < 4:
        # Need at least 4 taxa for 2 strata of size 2
        strata = {name: 0 for name in names}
        if not return_diagnostics:
            return strata
        return strata, {
            "n_strata": 1,
            "n_strata_natural": 1,
            "natural_bin_ids": np.zeros(n, dtype=np.intp),
            "min_stratum_size_used": min_stratum_size,
            "isolation": np.zeros(n, dtype=np.float64),
        }

    enc, is_res = _encode_msa(msa)
    D = _pairwise_hamming(enc, is_res)
    isolation = _min_nn_distances(D)

    # Cap strata at 2: captures the primary evolutionary split without
    # running into small-sample issues for permutation.
    max_clusters = min(2, max_clusters) if max_clusters else 2

    # Always compute natural strata first (unconstrained, min_s=2)
    natural_ids = _auto_strata(isolation, min_stratum_size=2, max_clusters=max_clusters)

    # If caller needs power-constrained strata, compute those too
    if min_stratum_size > 2:
        bin_ids = _auto_strata(isolation, min_stratum_size, max_clusters)
    else:
        bin_ids = natural_ids

    strata = {name: int(bin_ids[i]) for i, name in enumerate(names)}

    if not return_diagnostics:
        return strata

    return strata, {
        "n_strata": int(bin_ids.max()) + 1,
        "n_strata_natural": int(natural_ids.max()) + 1,
        "natural_bin_ids": natural_ids,
        "min_stratum_size_used": min_stratum_size,
        "isolation": isolation,
    }
