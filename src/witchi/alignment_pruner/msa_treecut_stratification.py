"""
MSA Isolation Stratification
=============================

Assigns aligned sequences to strata based on evolutionary isolation
(mean k-nearest-neighbour distance), for use in stratified permutation
testing of compositional bias.

Pipeline
--------
1. **Seed selection** — farthest-point sampling picks k well-spaced
   sequences (k ~ log2 N, typically 3-10).           O(k*N*L)
2. **Seed-distance projection** — gap-aware Hamming distance from each
   seed to every sequence; each sequence becomes a point in R^k.
   Streaming chunks avoid materialising the full N*L matrix.
                                                      O(k*N*L)
3. **Nearest-neighbour graph** — in the low-dimensional projection
   space, build a k-NN graph via dimension sweeps.    O(N*d*log N)
4. **Isolation metric** — mean k-NN distance per taxon, a proxy for
   terminal branch length / compositional drift.      O(N)
5. **Unsupervised binning** — Fisher-Jenks natural breaks on the 1-D
   isolation metric, with Calinski-Harabasz index for automatic K
   selection.  No thresholds.                         O(N*K_max^2)

Total: **O(k*N*L + N*log N)** — nearly linear in the MSA size N*L.

Designed for phylogenomic concatenated alignments (L up to 100,000+,
N up to thousands).

Requirements: numpy only.
"""

from __future__ import annotations

import math
from collections import defaultdict
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
# Phase 1-2: Seed projection  — O(k*N*L), memory O(N*L + N*k)
# ---------------------------------------------------------------------------

def _hamming_seed_to_all(enc: np.ndarray, is_res: np.ndarray,
                         seed_idx: int) -> np.ndarray:
    """Fractional Hamming from one seed to all sequences (fully vectorised).

    Parameters
    ----------
    enc : (N, L) uint8 — pre-encoded MSA
    is_res : (N, L) bool — residue mask (True = not a gap)
    seed_idx : index of seed sequence
    """
    shared = is_res & is_res[seed_idx]                      # (N, L) bool
    n_shared = shared.sum(axis=1).astype(np.float64)         # (N,)
    mismatches = ((enc != enc[seed_idx]) & shared).sum(      # (N,)
        axis=1
    ).astype(np.float64)
    return np.where(n_shared > 0, mismatches / n_shared, 1.0)


def _seed_projection(msa: list[str], k: Optional[int] = None,
                     rng_seed: int = 42) -> np.ndarray:
    """Project sequences into R^k via distances to k seed sequences.

    Encodes the MSA once, selects seeds via farthest-point sampling,
    and reuses seed distances computed during selection (no redundant
    second pass).  Returns shape (N, k).
    """
    n = len(msa)
    if k is None:
        k = max(3, min(int(math.ceil(math.log2(max(n, 4)))), 12))
    k = min(k, n)

    enc, is_res = _encode_msa(msa)
    rng = np.random.default_rng(rng_seed)

    # Farthest-point sampling — save per-seed distances for the projection
    seeds = [int(rng.integers(n))]
    seed_dists = [_hamming_seed_to_all(enc, is_res, seeds[0])]
    min_dist = seed_dists[0].copy()

    for _ in range(k - 1):
        cand = min_dist.copy()
        for s in seeds:
            cand[s] = -1.0
        new = int(np.argmax(cand))
        seeds.append(new)
        d = _hamming_seed_to_all(enc, is_res, new)
        seed_dists.append(d)
        min_dist = np.minimum(min_dist, d)

    return np.column_stack(seed_dists)  # (N, k)


# ---------------------------------------------------------------------------
# Phase 3: Nearest-neighbour graph  — O(N*d*log N)
# ---------------------------------------------------------------------------

def _euclidean(a: np.ndarray, b: np.ndarray) -> float:
    d = a - b
    return float(math.sqrt(d @ d))


def _build_nn_graph(proj: np.ndarray, k_nn: int = 6
                    ) -> list[tuple[int, int, float]]:
    """k-NN graph via dimension sweeps in low-d projection space."""
    n, d = proj.shape
    k_nn = min(k_nn, n - 1)

    if n <= 800:
        return _nn_graph_brute(proj, k_nn)

    edges: dict[tuple[int, int], float] = {}
    window = 2 * k_nn + 1

    for dim in range(d):
        order = np.argsort(proj[:, dim])
        for pos in range(n):
            i = int(order[pos])
            lo = max(0, pos - window)
            hi = min(n, pos + window + 1)
            for q in range(lo, hi):
                j = int(order[q])
                if i >= j:
                    continue
                key = (i, j)
                if key not in edges:
                    edges[key] = _euclidean(proj[i], proj[j])

    return _keep_k_nearest(edges, k_nn)


def _nn_graph_brute(proj: np.ndarray, k_nn: int
                    ) -> list[tuple[int, int, float]]:
    """Exact k-NN for small N — vectorised pairwise distances."""
    n = proj.shape[0]
    diff = proj[:, np.newaxis, :] - proj[np.newaxis, :, :]
    dists = np.sqrt((diff * diff).sum(axis=2))

    edges: dict[tuple[int, int], float] = {}
    for i in range(n):
        order = np.argsort(dists[i])
        for j in order[1:k_nn + 1]:
            key = (min(i, int(j)), max(i, int(j)))
            edges[key] = float(dists[i, j])
    return [(u, v, w) for (u, v), w in edges.items()]


def _keep_k_nearest(edges: dict[tuple[int, int], float], k: int
                    ) -> list[tuple[int, int, float]]:
    per_node: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for (u, v), w in edges.items():
        per_node[u].append((w, v))
        per_node[v].append((w, u))
    kept: dict[tuple[int, int], float] = {}
    for node, nbrs in per_node.items():
        nbrs.sort()
        for w, other in nbrs[:k]:
            key = (min(node, other), max(node, other))
            kept[key] = w
    return [(u, v, w) for (u, v), w in kept.items()]


# ---------------------------------------------------------------------------
# Phase 4: Isolation metric  — O(N)
# ---------------------------------------------------------------------------

def _mean_knn_distances(
    n: int,
    edges: list[tuple[int, int, float]],
) -> np.ndarray:
    """Mean k-NN distance per taxon from the k-NN graph.

    Returns shape (N,) array of mean distances.
    """
    per_node: dict[int, list[float]] = defaultdict(list)
    for u, v, w in edges:
        per_node[u].append(w)
        per_node[v].append(w)
    out = np.zeros(n, dtype=np.float64)
    for i in range(n):
        dists = per_node.get(i, [])
        out[i] = float(np.mean(dists)) if dists else 0.0
    return out


# ---------------------------------------------------------------------------
# Phase 5: Fisher-Jenks natural breaks + Calinski-Harabasz  — O(N*K^2)
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
    n_neighbours: int = 6,
    n_seeds: Optional[int] = None,
    return_diagnostics: bool = False,
) -> dict[str, int] | tuple[dict[str, int], dict]:
    """Assign each sequence to a stratum for stratified permutation testing.

    Strata group taxa by evolutionary isolation (mean k-NN distance),
    so that within-stratum permutations preserve the correlation structure
    imposed by shared evolutionary history.  The number of strata is
    selected automatically via Calinski-Harabasz on Fisher-Jenks natural
    breaks.

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
    n_neighbours : int
        NN graph connectivity (default 6).
    n_seeds : int, optional
        Seed count for distance projection.
        Default: clamp(ceil(log2 N), 3, 12).
    return_diagnostics : bool
        If True, return (strata_dict, diagnostics_dict) instead of just
        strata_dict.  Diagnostics include the isolation metric and the
        number of natural strata (min_stratum_size=2) for comparison.

    Returns
    -------
    dict mapping sequence name (str) -> stratum index (int)
        — or (dict, diagnostics) if return_diagnostics=True
    """
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

    proj = _seed_projection(msa, k=n_seeds)
    edges = _build_nn_graph(proj, k_nn=n_neighbours)
    isolation = _mean_knn_distances(n, edges)

    # Always compute natural strata first (unconstrained, min_s=2, K_max capped at 10)
    natural_max_k = min(10, max_clusters) if max_clusters else 10
    natural_ids = _auto_strata(isolation, min_stratum_size=2, max_clusters=natural_max_k)

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
