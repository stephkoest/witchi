"""
Similarity-stratified permutation testing for WitChi.

Assigns taxa to strata based on evolutionary isolation (min nearest-neighbour
distance under gap-aware Hamming, a direct proxy for terminal branch length)
and performs global-baseline column permutation within strata.  This preserves
the correlation structure imposed by shared evolutionary history on uneven
trees, producing a null distribution that accounts for branch-length-driven
compositional drift.

Standard permutation is the special case where all taxa belong to one
stratum.

Key design choices
------------------
- **Per-stratum chi2**: each taxon's chi2 is computed against its own
  stratum's expected frequencies, not global.  This eliminates the
  stratum-global offset that causes Jensen's-inequality inflation when
  combined with within-stratum permutation (observed and null chi2 share
  the same reference frame).
- **Per-stratum p-values**: each taxon's p-value is computed against its
  own stratum's null pool, then Bonferroni-corrected by N.  Prevents
  cross-stratum contamination from fat-tailed long-branch strata.
- **Mean-centred robust Z-scores**: Z = (x - mean) / (1.4826 * MAD).
  Mean centering is required because the mean has the tower property
  (consistent across conditional and marginal nulls) while the median
  does not.
- **Global chi2 for pruning**: column scoring (delta-chi2) still uses
  global expected frequencies — pruning needs a common reference across
  all taxa.
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

    def __init__(
        self,
        sums,
        maxes,
        upper_box_threshold,
        upper_threshold,
        pooled_null,
        chi2_matrix,
        stratum_pools,
        bin_ids,
        diagnostics,
    ):
        # Standard fields (same semantics as PermutationTest.run() output)
        self.sums = sums  # (P,) alignment-level chi2 sums
        self.maxes = maxes  # (P,) per-permutation max chi2
        self.upper_box_threshold = upper_box_threshold  # 75th percentile of pooled null
        self.upper_threshold = upper_threshold  # 95th percentile of pooled null
        self.pooled_null = pooled_null  # (P*N,) flattened pooled null

        # Stratified-specific fields
        self.chi2_matrix = chi2_matrix  # (P, N) per-taxon chi2 per permutation
        self.stratum_pools = stratum_pools  # {taxon_idx: pool_array}
        self.bin_ids = bin_ids  # (N,) stratum assignment
        self.diagnostics = diagnostics  # natural/realizable strata info

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
# Main entry point
# ---------------------------------------------------------------------------


def run_similarity_stratified(
    alignment_array,
    alignment,
    chi_square_calculator,
    permutation_test,
    max_clusters=None,
    n_neighbours=None,  # deprecated; isolation is exact min-NN
    n_seeds=None,  # deprecated; no projection step
):
    """Run similarity-stratified permutation test.

    Assigns taxa to strata by evolutionary isolation (min nearest-neighbour
    distance under gap-aware Hamming), then performs global-baseline column
    permutation within strata.

    Parameters
    ----------
    alignment_array : (N, L) numpy char array
        The alignment as read by AlignmentReader.
    alignment : Bio.Align.MultipleSeqAlignment
        BioPython alignment object (used for taxon names and raw sequences).
    chi_square_calculator : ChiSquareCalculator
        Configured chi2 calculator for this alignment's character set.
    permutation_test : PermutationTest
        Instance used for the permutation loop and summary statistics.
    max_clusters : int or None
        Upper bound on strata count.  Default: min(20, N // min_stratum_size).
    n_neighbours : int, optional
        Deprecated.  Kept for API compatibility; ignored.
    n_seeds : int, optional
        Deprecated.  Kept for API compatibility; ignored.

    Returns
    -------
    StratifiedResult
        Contains the standard 5-tuple fields (via .as_standard_tuple()) plus
        chi2_matrix, stratum_pools, bin_ids, and diagnostics for per-stratum
        p-value computation and reporting.
    """
    N = alignment_array.shape[0]
    num_permutations = permutation_test.permutations
    taxon_names = [rec.id for rec in alignment]
    raw_seqs = [str(rec.seq) for rec in alignment]

    # --- Derive min stratum size from N and P for Bonferroni power ---
    # Need: N / (s * P) < 0.05  =>  s > 20*N/P
    auto_min = max(2, math.ceil(20 * N / num_permutations) + 1)
    print(
        f"[similarity_stratified] auto min_stratum_size={auto_min} "
        f"(N={N}, P={num_permutations})"
    )

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
        p_recommended = math.ceil(20 * N / (s_min - 1))

    # Report strata
    print(
        f"[similarity_stratified] natural strata: {n_natural}, "
        f"realizable strata (at P={num_permutations}): {n_realizable}"
    )
    if n_natural > n_realizable and p_recommended is not None:
        print(
            f"[similarity_stratified] NOTE: {n_natural} natural strata "
            f"detected but only {n_realizable} realizable at "
            f"P={num_permutations}. To realize all strata, use "
            f"P>={p_recommended}."
        )

    # Convert name->stratum dict to index arrays
    n_strata = max(name_to_stratum.values()) + 1
    bin_ids = np.array([name_to_stratum[n] for n in taxon_names], dtype=int)
    bins = [np.where(bin_ids == k)[0] for k in range(n_strata)]
    bins = [b for b in bins if len(b) > 0]
    bin_sizes = [len(b) for b in bins]

    print(
        f"[similarity_stratified] {len(bins)} realizable strata, " f"sizes={bin_sizes}"
    )

    # --- Phase 2: global-baseline stratified permutation ---
    chi2_matrix = permutation_test._permute_and_calculate_chi2(
        alignment_array,
        chi_square_calculator,
        strata_indices=bins,
    )

    # --- Phase 3: build per-stratum null pools ---
    stratum_pools = {}
    for bin_indices in bins:
        pool = chi2_matrix[:, bin_indices].ravel().astype(np.float64)
        for idx in bin_indices:
            stratum_pools[int(idx)] = pool

    # --- Phase 4: standard summary statistics ---
    sums, maxes, upper_box_threshold, upper_threshold, pooled_null = (
        permutation_test._summarize_null(chi2_matrix)
    )

    diagnostics = {
        "n_strata_natural": n_natural,
        "n_strata_realizable": n_realizable,
        "bin_sizes": bin_sizes,
        "natural_stratum_sizes": natural_stratum_sizes,
        "p_recommended": p_recommended,
        "min_stratum_size_auto": auto_min,
        "name_to_stratum": {str(k): int(v) for k, v in name_to_stratum.items()},
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
