# stratified_permutation.py
import numpy as np
from .permutation_test import PermutationTest
from .chi_square_calculator import ChiSquareCalculator


def calc_per_row_chi2(alignment_array: np.ndarray, chi_calc: ChiSquareCalculator):
    """
    Reproduce your per-row χ² computation using WitChi calculator:
      counts -> expected -> per-row χ²; plus alignment sum.
    """
    count_rows = chi_calc.calculate_row_counts(alignment_array)
    expected = chi_calc.calculate_expected_observed(count_rows)
    per_row = chi_calc.calculate_row_chi2(expected, count_rows)
    return per_row, float(np.sum(per_row))


def _adaptive_stratified_bins(
    obs_chi: np.ndarray,
    target_K: int = 4,
    K_min: int = 2,
    K_max: int = 20,
    min_bin_frac: float = 0.01,
):
    """Return quantile bin edges that avoid tiny bins.

    - Start from target_K; if any bin would contain < min_bin (as a fraction of N),
      shrink K until all bins are sufficiently populated (but not below K_min).
    - If all bins are very large and flat, allow up to K_max by splitting until
      a bin would violate the minimum.

    Returns:
      qs: 1D array of edges (length K+1)
      K:  effective number of bins
      counts: bin counts for the final partition (length K)
    """
    N = int(obs_chi.size)
    if N == 0:
        return np.array([0.0, 1.0]), 1, np.array([0], dtype=int)
    min_bin = max(2, int(np.ceil(min_bin_frac * N)))

    def edges_and_counts(K):
        qs = np.quantile(obs_chi, np.linspace(0, 1, K + 1))
        ids = np.clip(np.searchsorted(qs, obs_chi, side="right") - 1, 0, K - 1)
        counts = np.bincount(ids, minlength=K)
        return qs, counts

    K = int(np.clip(target_K, K_min, K_max))
    qs, counts = edges_and_counts(K)

    # shrink if any tiny bins
    while K > K_min and np.any(counts < min_bin):
        K -= 1
        qs, counts = edges_and_counts(K)

    # optionally expand a bit if comfortably above threshold
    while K < K_max:
        qs_next, counts_next = edges_and_counts(K + 1)
        if np.any(counts_next < min_bin):
            break
        K += 1
        qs, counts = qs_next, counts_next

    # concise console summary
    props = (counts / counts.sum()) if counts.sum() > 0 else counts.astype(float)
    print(
        f"[Stratum bins] K={K}  counts={counts.tolist()}  "
        f"props={[round(float(p),3) for p in props]}  "
        f"edges={qs.tolist()}"
    )
    return qs, K, counts


def run_stratified(
    alignment_array: np.ndarray,
    chi_calc: ChiSquareCalculator,
    permutations: int,
    K_bins: int = 4,
    adaptive_bins: bool = True,
    min_bin_frac: float = 0.01,
):
    """
    χ² stratified permutation (stratified): bin taxa by observed χ² (TPSB-style), then for each stratum,
    run column permutations ONLY within that stratum and pool all per-taxon χ².
    The total number of permutations is split across strata proportionally to stratum size.
    """
    obs_chi, _ = calc_per_row_chi2(alignment_array, chi_calc)
    if adaptive_bins:
        qs, K_use, init_counts = _adaptive_stratified_bins(
            obs_chi, target_K=K_bins, K_min=2, K_max=20, min_bin_frac=min_bin_frac
        )
    else:
        K_use = K_bins
        qs = np.quantile(obs_chi, np.linspace(0, 1, K_use + 1))
    bin_ids = np.clip(np.searchsorted(qs, obs_chi, side="right") - 1, 0, K_use - 1)
    bins = [np.where(bin_ids == k)[0] for k in range(K_use)]
    bins = [b for b in bins if len(b) > 0]
    bin_sizes = np.array([len(b) for b in bins], dtype=int)

    pooled = []
    skipped_singleton_bins = 0
    perm_engine_bin = PermutationTest(
        num_workers_permute=1, permutations=int(permutations)
    )
    for b_idx, (bin_indices, sz) in enumerate(zip(bins, bin_sizes)):
        if sz < 2:
            skipped_singleton_bins += 1
            continue
        subalignment = alignment_array[bin_indices, :]
        chi_mat = perm_engine_bin._permute_and_calculate_chi2(subalignment, chi_calc)
        # chi_mat shape: (perm_count, sz)
        pooled.append(np.asarray(chi_mat, dtype=np.float64).ravel())

    # pooled_arr = np.array(pooled, dtype=np.float64)
    pooled_arr = np.concatenate(pooled) if pooled else np.array([], dtype=np.float64)
    return pooled_arr, {"K": K_use, "sizes": bin_sizes}
