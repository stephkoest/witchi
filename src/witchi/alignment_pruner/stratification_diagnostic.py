"""Diagnostic for similarity-stratified permutation validity.

Compares alignment-level chi2 sum distributions from standard vs
stratified permutation using three metrics:

1. **Concordance** (hard gate, threshold-free): do the observed
   alignment chi2's empirical p-values under both nulls lead to the
   same significance conclusion?  If the standard null says biased but
   the stratified null does not, the inflation is masking real signal.

2. **Inflation Z** (robust Z-score): how far the stratified null median
   is from the standard null, in robust standard deviations (MAD/0.6745).
   Z >= 3 is the conventional outlier threshold.

3. **Signal consumed** (relative inflation): what fraction of the
   observed signal is consumed by the null inflation.  Puts the absolute
   inflation into perspective — a large Z is harmless if the observed
   chi2 dwarfs both nulls.
"""

import numpy as np

from .utils import _robust_zscore


_ZSCORE_OUTLIER_THRESHOLD = 3.0


# ---------------------------------------------------------------------------
# Core comparison (reusable with pre-computed data)
# ---------------------------------------------------------------------------

def compare_null_distributions(observed_chi2, std_sums, strat_sums,
                               n_strata, alpha=0.05):
    """Compare standard vs stratified null chi2 sum distributions.

    Parameters
    ----------
    observed_chi2 : float
        Observed alignment-level chi2 sum.
    std_sums : (P,) array
        Alignment-level chi2 sums from standard permutation.
    strat_sums : (P,) array
        Alignment-level chi2 sums from stratified permutation.
    n_strata : int
        Number of realizable strata.
    alpha : float
        Significance level for concordance check.

    Returns
    -------
    dict
        ``valid``              – True if p-values are concordant.
        ``warning``            – True if inflation_z >= 3 AND
                                 signal_consumed >= 0.10.
        ``inflation_z``        – robust Z-score of the stratified null
                                 median against the standard null.
        ``signal_consumed``    – fraction of observed signal consumed by
                                 the inflation: (strat_med - std_med) /
                                 (observed - std_med).  0 = no impact,
                                 1 = inflation equals the entire signal.
        ``observed_chi2``      – echo of input.
        ``p_standard``         – empirical p under standard null.
        ``p_stratified``       – empirical p under stratified null.
        ``concordant``         – same as ``valid``.
        ``standard_median``    – median of standard sums.
        ``stratified_median``  – median of stratified sums.
        ``n_strata``           – echo of input.
        ``alpha``              – echo of input.
    """
    p_std = float(np.sum(std_sums >= observed_chi2) / len(std_sums))
    p_strat = float(np.sum(strat_sums >= observed_chi2) / len(strat_sums))

    concordant = (p_std < alpha) == (p_strat < alpha)

    std_med = float(np.median(std_sums))
    strat_med = float(np.median(strat_sums))
    inflation_z = float(_robust_zscore(strat_med, std_sums))

    # Fraction of observed signal consumed by the inflation
    signal_gap = observed_chi2 - std_med
    if signal_gap > 0:
        signal_consumed = (strat_med - std_med) / signal_gap
    else:
        # Observed <= standard null median: no signal to consume
        signal_consumed = 0.0

    # Warning requires both: statistically outlying AND practically relevant
    warning = inflation_z >= _ZSCORE_OUTLIER_THRESHOLD and signal_consumed >= 0.10

    return {
        "valid": concordant,
        "warning": warning,
        "inflation_z": inflation_z,
        "signal_consumed": float(signal_consumed),
        "observed_chi2": float(observed_chi2),
        "p_standard": p_std,
        "p_stratified": p_strat,
        "concordant": concordant,
        "standard_median": std_med,
        "stratified_median": strat_med,
        "n_strata": n_strata,
        "alpha": alpha,
    }


def print_diagnostic(diag):
    """Print a human-readable summary of the diagnostic dict."""
    print(
        f"[diagnostic] Observed chi2: {diag['observed_chi2']:.2f} | "
        f"P(standard): {diag['p_standard']:.4f} | "
        f"P(stratified): {diag['p_stratified']:.4f}"
    )
    print(
        f"[diagnostic] Standard null median: {diag['standard_median']:.2f} | "
        f"Stratified null median: {diag['stratified_median']:.2f} | "
        f"Inflation Z: {diag['inflation_z']:.1f} | "
        f"Signal consumed: {diag['signal_consumed']:.1%}"
    )
    if not diag["valid"]:
        print(
            "[diagnostic] WARNING: p-values are discordant — stratified "
            "null inflation changes the significance call. Stratification "
            "is not valid for this alignment."
        )
    elif diag["warning"]:
        print(
            "[diagnostic] CAUTION: null distributions are concordant but "
            "inflation consumes >= 10% of the observed signal. Z-scores "
            "and pruning depth may be affected."
        )
    else:
        print("[diagnostic] Stratification valid — null distributions concordant.")


# ---------------------------------------------------------------------------
# Standalone entry point (runs both permutations from scratch)
# ---------------------------------------------------------------------------

def diagnose_stratification_validity(
    alignment_array,
    alignment,
    chi_square_calculator,
    permutations=100,
    num_workers=1,
    alpha=0.05,
    max_clusters=None,
    n_neighbours=6,
    n_seeds=None,
):
    """Diagnose whether similarity-stratified permutation inflates the null.

    Runs both standard and stratified permutation on the same alignment,
    computes the observed alignment chi2 sum, and checks whether its
    empirical p-value under both nulls leads to the same conclusion at
    the given *alpha*.

    Parameters
    ----------
    alignment_array : (N, L) numpy char array
    alignment : Bio.Align.MultipleSeqAlignment
    chi_square_calculator : ChiSquareCalculator
    permutations : int
        Number of permutations for both runs.
    num_workers : int
        Parallel workers for permutation.
    alpha : float
        Significance level used for concordance check.
    max_clusters, n_neighbours, n_seeds :
        Passed through to ``run_similarity_stratified``.

    Returns
    -------
    dict
        See :func:`compare_null_distributions` for keys.  Adds
        ``trivial`` (bool) — True if only 1 stratum.
    """
    from .permutation_test import PermutationTest
    from .stratified_permutation import run_similarity_stratified

    # --- Observed alignment chi2 ---
    row_counts = chi_square_calculator.calculate_row_counts(alignment_array)
    exp_obs = chi_square_calculator.calculate_expected_observed(row_counts)
    per_row_chi2 = chi_square_calculator.calculate_row_chi2(exp_obs, row_counts)
    observed_chi2 = float(np.sum(per_row_chi2))

    # --- Stratified permutation (also reveals n_strata) ---
    pt_strat = PermutationTest(
        num_workers_permute=num_workers, permutations=permutations,
    )
    strat_result = run_similarity_stratified(
        alignment_array, alignment, chi_square_calculator,
        permutation_test=pt_strat,
        max_clusters=max_clusters,
        n_neighbours=n_neighbours,
        n_seeds=n_seeds,
    )
    n_strata = strat_result.diagnostics["n_strata_realizable"]
    strat_sums = strat_result.sums

    # --- Trivial case (1 stratum = identical to standard) ---
    if n_strata <= 1:
        p = float(np.sum(strat_sums >= observed_chi2) / len(strat_sums))
        med = float(np.median(strat_sums))
        result = {
            "valid": True,
            "warning": False,
            "inflation_z": 0.0,
            "signal_consumed": 0.0,
            "observed_chi2": observed_chi2,
            "p_standard": p,
            "p_stratified": p,
            "concordant": True,
            "standard_median": med,
            "stratified_median": med,
            "n_strata": 1,
            "trivial": True,
            "alpha": alpha,
        }
        return result

    # --- Standard permutation ---
    pt_std = PermutationTest(
        num_workers_permute=num_workers, permutations=permutations,
    )
    std_sums, _, _, _, _ = pt_std.run(alignment_array, chi_square_calculator)

    result = compare_null_distributions(
        observed_chi2, std_sums, strat_sums, n_strata, alpha=alpha,
    )
    result["trivial"] = False
    return result
