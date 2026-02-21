"""
Similarity-stratified permutation testing for WitChi.

Assigns taxa to strata based on evolutionary isolation (mean k-NN distance
in a seed-projection space) and performs global-baseline column permutation
within strata.  This preserves the correlation structure imposed by shared
evolutionary history on uneven trees, producing a null distribution that
accounts for branch-length-driven compositional drift.

Standard permutation is the special case where all taxa belong to one
stratum.

Key design choices
------------------
- **Global-baseline permutation**: permute within strata on the full
  alignment, compute chi2 against global expected frequencies.  Observed
  and null chi2 measure the same thing (deviation from global composition).
- **Per-stratum p-values**: each taxon's p-value is computed against its
  own stratum's null pool, then Bonferroni-corrected by N.  Prevents
  cross-stratum contamination from fat-tailed long-branch strata.
- **Pooled null for pruning**: the full (P*N) pooled null is used as the
  Wasserstein pruning target.  Its fat-tailed shape correctly represents
  the true null on an uneven tree.
- **Natural vs realizable strata**: two clusterings are computed.  Natural
  strata (min_s=2) reveal the true isolation structure.  Realizable strata
  (min_s = ceil(20*N/P)+1) ensure Bonferroni-corrected significance is
  achievable.  Realizable strata are used for permutation.
"""

import math

import numpy as np
from joblib import Parallel, delayed

from .chi_square_calculator import ChiSquareCalculator
from .msa_treecut_stratification import msa_strata


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

class StratifiedResult:
    """Result from similarity-stratified permutation testing.

    Holds the standard fields produced by PermutationTest.run() (sums,
    maxes, upper_box_threshold, upper_threshold, pooled_null) so the
    pruning loop can consume them unchanged, plus stratified-specific
    fields for per-stratum p-value computation.
    """

    def __init__(self, sums, maxes, upper_box_threshold, upper_threshold,
                 pooled_null, chi2_matrix, stratum_pools, bin_ids, diagnostics):
        # Standard fields (same semantics as PermutationTest.run() output)
        self.sums = sums                          # (P,) alignment-level chi2 sums
        self.maxes = maxes                        # (P,) per-permutation max chi2
        self.upper_box_threshold = upper_box_threshold  # 75th percentile of pooled null
        self.upper_threshold = upper_threshold    # 95th percentile of pooled null
        self.pooled_null = pooled_null            # (P*N,) flattened pooled null

        # Stratified-specific fields
        self.chi2_matrix = chi2_matrix            # (P, N) per-taxon chi2 per permutation
        self.stratum_pools = stratum_pools        # {taxon_idx: pool_array}
        self.bin_ids = bin_ids                    # (N,) stratum assignment
        self.diagnostics = diagnostics            # natural/realizable strata info

    def as_standard_tuple(self):
        """Unpack to the 5-tuple matching PermutationTest.run() output.

        Returns (sums, maxes, upper_box_threshold, upper_threshold, pooled_null).
        """
        return (
            self.sums,
            self.maxes,
            self.upper_box_threshold,
            self.upper_threshold,
            self.pooled_null,
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def calc_per_row_chi2(alignment_array, chi_calc):
    """Compute per-taxon chi2 and alignment sum using WitChi calculator.

    Returns (per_row_chi2, alignment_sum).
    """
    count_rows = chi_calc.calculate_row_counts(alignment_array)
    expected = chi_calc.calculate_expected_observed(count_rows)
    per_row = chi_calc.calculate_row_chi2(expected, count_rows)
    return per_row, float(np.sum(per_row))


def _permute_within_strata(alignment_array, strata_indices, chi_calc,
                           permutations, num_workers=1):
    """Global-baseline stratified column permutation.

    For each permutation replicate:
      1. Copy the full N x L alignment.
      2. For each stratum, independently shuffle each column's residues
         among that stratum's taxa.
      3. Compute per-taxon chi2 on the full permuted alignment using
         global expected frequencies.

    Parameters
    ----------
    alignment_array : (N, L) char array
    strata_indices : list of 1-D int arrays, one per stratum
    chi_calc : ChiSquareCalculator
    permutations : int
    num_workers : int
        Joblib parallel workers (matches PermutationTest convention).

    Returns
    -------
    chi2_matrix : (P, N) float64 array
    alignment_sums : (P,) float64 array
    """
    def _single_permutation(i):
        rng = np.random.default_rng(12345 + i)
        permuted = alignment_array.copy()
        for s_indices in strata_indices:
            if len(s_indices) < 2:
                continue
            sub = permuted[s_indices, :]
            permuted[s_indices, :] = np.apply_along_axis(
                rng.permutation, 0, sub
            )
        count_rows = chi_calc.calculate_row_counts(permuted)
        expected = chi_calc.calculate_expected_observed(count_rows)
        per_row = chi_calc.calculate_row_chi2(expected, count_rows)
        return per_row

    results = Parallel(n_jobs=num_workers)(
        delayed(_single_permutation)(i) for i in range(permutations)
    )

    chi2_matrix = np.array(results, dtype=np.float64)
    alignment_sums = np.sum(chi2_matrix, axis=1)
    return chi2_matrix, alignment_sums


# ---------------------------------------------------------------------------
# Per-stratum p-value computation
# ---------------------------------------------------------------------------

def calc_empirical_pvalue_per_stratum(per_row_chi2, stratum_pools, n_taxa):
    """Compute Bonferroni-corrected empirical p-values using per-stratum pools.

    Each taxon's p-value is the fraction of its stratum's null pool >= the
    observed chi2, multiplied by n_taxa (Bonferroni correction), clipped
    to [0, 1].

    Parameters
    ----------
    per_row_chi2 : (N,) array of observed per-taxon chi2
    stratum_pools : dict mapping taxon index (int) -> (s*P,) null pool array
    n_taxa : int, total number of taxa (N) for Bonferroni correction

    Returns
    -------
    list of float, length N
    """
    pvals = []
    for i in range(len(per_row_chi2)):
        pool = stratum_pools[i]
        p = float(np.sum(pool >= per_row_chi2[i]) / len(pool)) * n_taxa
        pvals.append(min(max(p, 0.0), 1.0))
    return pvals


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_similarity_stratified(
    alignment_array,
    alignment,
    chi_square_calculator,
    num_permutations,
    num_workers=1,
    max_clusters=None,
    n_neighbours=6,
    n_seeds=None,
):
    """Run similarity-stratified permutation test.

    Assigns taxa to strata by evolutionary isolation (mean k-NN distance),
    then performs global-baseline column permutation within strata.

    Parameters
    ----------
    alignment_array : (N, L) numpy char array
        The alignment as read by AlignmentReader.
    alignment : Bio.Align.MultipleSeqAlignment
        BioPython alignment object (used for taxon names and raw sequences).
    chi_square_calculator : ChiSquareCalculator
        Configured chi2 calculator for this alignment's character set.
    num_permutations : int
        Number of permutation replicates (P).
    num_workers : int
        Parallel workers for the permutation loop.
    max_clusters : int or None
        Upper bound on strata count.  Default: min(20, N // min_stratum_size).
    n_neighbours : int
        k-NN connectivity for the isolation metric (default 6).
    n_seeds : int or None
        Seed count for distance projection.  Default: auto.

    Returns
    -------
    StratifiedResult
        Contains the standard 5-tuple fields (via .as_standard_tuple()) plus
        chi2_matrix, stratum_pools, bin_ids, and diagnostics for per-stratum
        p-value computation and reporting.
    """
    N = alignment_array.shape[0]
    taxon_names = [rec.id for rec in alignment]
    raw_seqs = [str(rec.seq) for rec in alignment]

    # --- Derive min stratum size from N and P for Bonferroni power ---
    # Need: N / (s * P) < 0.05  =>  s > 20*N/P
    auto_min = max(2, math.ceil(20 * N / num_permutations) + 1)
    print(f"[similarity_stratified] auto min_stratum_size={auto_min} "
          f"(N={N}, P={num_permutations})")

    # --- Phase 1: compute strata from evolutionary isolation ---
    name_to_stratum, diag = msa_strata(
        msa=raw_seqs,
        names=taxon_names,
        min_stratum_size=auto_min,
        max_clusters=max_clusters,
        n_neighbours=n_neighbours,
        n_seeds=n_seeds,
        return_diagnostics=True,
    )

    n_natural = diag["n_strata_natural"]
    n_realizable = diag["n_strata"]
    natural_ids = diag["natural_bin_ids"]

    # P recommendation from natural strata
    p_recommended = None
    natural_stratum_sizes = None
    if n_natural > 1:
        natural_stratum_sizes = [
            int((natural_ids == k).sum()) for k in range(n_natural)
        ]
        s_min = min(natural_stratum_sizes)
        p_recommended = math.ceil(20 * N / s_min) + 1

    # Report strata
    print(f"[similarity_stratified] natural strata: {n_natural}, "
          f"realizable strata (at P={num_permutations}): {n_realizable}")
    if n_natural > n_realizable and p_recommended is not None:
        print(f"[similarity_stratified] NOTE: {n_natural} natural strata "
              f"detected but only {n_realizable} realizable at "
              f"P={num_permutations}. To realize all strata, use "
              f"P>={p_recommended}.")

    # Convert name->stratum dict to index arrays
    n_strata = max(name_to_stratum.values()) + 1
    bin_ids = np.array(
        [name_to_stratum[n] for n in taxon_names], dtype=int
    )
    bins = [np.where(bin_ids == k)[0] for k in range(n_strata)]
    bins = [b for b in bins if len(b) > 0]
    bin_sizes = [len(b) for b in bins]

    print(f"[similarity_stratified] {len(bins)} realizable strata, "
          f"sizes={bin_sizes}")
    print(f"Running {num_permutations} permutations.")
    print(f"Using {num_workers} worker(s) for permutation")

    # --- Phase 2: global-baseline stratified permutation ---
    chi2_matrix, alignment_sums = _permute_within_strata(
        alignment_array, bins, chi_square_calculator,
        num_permutations, num_workers,
    )

    # --- Phase 3: build per-stratum null pools ---
    # Map each taxon index to the pool of its stratum's chi2 values.
    # Pool for stratum s = chi2_matrix[:, stratum_indices].ravel()
    #   → s_i * P values, naturally weighted by stratum size.
    stratum_pools = {}
    for bin_indices in bins:
        pool = chi2_matrix[:, bin_indices].ravel().astype(np.float64)
        for idx in bin_indices:
            stratum_pools[int(idx)] = pool

    # --- Phase 4: standard summary statistics ---
    # Match PermutationTest.run() output format so the pruning loop
    # can consume the pooled null unchanged.
    pooled_null = chi2_matrix.ravel().astype(np.float64)
    sums = alignment_sums
    maxes = np.max(chi2_matrix, axis=1)
    upper_box_threshold = float(np.percentile(pooled_null, 75))
    upper_threshold = float(np.percentile(pooled_null, 95))

    diagnostics = {
        "n_strata_natural": n_natural,
        "n_strata_realizable": n_realizable,
        "bin_sizes": bin_sizes,
        "natural_stratum_sizes": natural_stratum_sizes,
        "p_recommended": p_recommended,
        "min_stratum_size_auto": auto_min,
        "name_to_stratum": {
            str(k): int(v) for k, v in name_to_stratum.items()
        },
    }

    return StratifiedResult(
        sums=sums,
        maxes=maxes,
        upper_box_threshold=upper_box_threshold,
        upper_threshold=upper_threshold,
        pooled_null=pooled_null,
        chi2_matrix=chi2_matrix,
        stratum_pools=stratum_pools,
        bin_ids=bin_ids,
        diagnostics=diagnostics,
    )
