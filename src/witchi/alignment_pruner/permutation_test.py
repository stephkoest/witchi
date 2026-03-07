import numpy as np
from joblib import Parallel, delayed

from .utils import write_score_dict_to_tsv, make_score_dict


def _freedman_lane_chi2(count_rows_array, strata_indices):
    """Compute per-taxon chi2 with Freedman-Lane recentering.

    Uses (observed - stratum_expected)^2 / global_expected per stratum.
    All strata are guaranteed to have >= 2 taxa (enforced by auto_min).
    """
    C, N = count_rows_array.shape
    total_per_char = count_rows_array.sum(axis=1)
    grand_total = total_per_char.sum()
    global_freq = total_per_char / grand_total
    taxon_totals = count_rows_array.sum(axis=0)
    global_expected = global_freq[:, np.newaxis] * taxon_totals[np.newaxis, :]
    safe_denom = np.where(global_expected > 0, global_expected, 1.0)

    per_row_chi2 = np.zeros(N, dtype=np.float64)

    for s_indices in strata_indices:
        stratum_counts = count_rows_array[:, s_indices]
        stratum_total = stratum_counts.sum(axis=1)
        stratum_freq = stratum_total / stratum_total.sum()
        stratum_expected = (
            stratum_freq[:, np.newaxis] * taxon_totals[s_indices][np.newaxis, :]
        )
        fl_chi2 = (stratum_counts - stratum_expected) ** 2 / safe_denom[:, s_indices]
        per_row_chi2[s_indices] = fl_chi2.sum(axis=0)

    return per_row_chi2


class PermutationTest:
    def __init__(self, num_workers_permute, permutations):
        self.is_dna = None
        self.chi_square_calculator = None
        self.num_workers = num_workers_permute
        self.permutations = permutations
        self._stratified_result = None

    def calc_empirical_pvalue(self, per_row_chi2, permutated_per_row_chi2,
                              per_taxon_pools=None):
        """Calculate empirical p-value.

        Parameters
        ----------
        per_row_chi2 : np.ndarray or scalar
            Observed per-taxon chi² (array) or alignment-level chi² (scalar).
        permutated_per_row_chi2 : np.ndarray
            Pooled null distribution (used when per_taxon_pools is None).
        per_taxon_pools : dict or None
            If provided, maps taxon index (int) to its stratum's null pool
            array.  Each taxon is compared against its own pool instead of
            the single pooled null.
        """
        empirical_p_list = []
        # check if array, indicating per taxon chi2 scores
        if isinstance(per_row_chi2, np.ndarray):
            for i in range(len(per_row_chi2)):
                pool = (per_taxon_pools[i] if per_taxon_pools is not None
                        else permutated_per_row_chi2)
                empirical_p = (
                    np.sum(pool >= per_row_chi2[i])
                    / len(pool)
                ) * len(per_row_chi2)
                empirical_p = min(max(empirical_p, 0.0), 1.0)
                empirical_p_list.append(empirical_p)
        # if not array, indicating total alignment chi2 scores
        else:
            empirical_p = np.sum(permutated_per_row_chi2 >= per_row_chi2) / len(
                permutated_per_row_chi2
            )
            empirical_p_list.append(empirical_p)

        return empirical_p_list

    def _permute_and_calculate_chi2(self, alignment_array, chi_square_calculator,
                                     strata_indices=None):
        def permute_single_calculate_chi2(i):
            iter_seed = 12345 + i
            rng = np.random.default_rng(iter_seed)
            if strata_indices is None:
                permuted_array = np.apply_along_axis(
                    rng.permutation, 0, alignment_array
                )
            else:
                permuted_array = alignment_array.copy()
                for s_indices in strata_indices:
                    if len(s_indices) < 2:
                        continue
                    sub = permuted_array[s_indices, :]
                    permuted_array[s_indices, :] = np.apply_along_axis(
                        rng.permutation, 0, sub
                    )
            count_rows_array = chi_square_calculator.calculate_row_counts(
                permuted_array
            )
            if strata_indices is not None:
                per_row_chi2 = _freedman_lane_chi2(
                    count_rows_array, strata_indices
                )
            else:
                expected_observed = (
                    chi_square_calculator.calculate_expected_observed(
                        count_rows_array
                    )
                )
                per_row_chi2 = chi_square_calculator.calculate_row_chi2(
                    expected_observed, count_rows_array
                )
            return per_row_chi2

        results = Parallel(n_jobs=self.num_workers)(
            delayed(permute_single_calculate_chi2)(i) for i in range(self.permutations)
        )
        return np.array(results, dtype=np.float64)

    @staticmethod
    def _summarize_null(chi2_matrix):
        """Compute summary statistics from a (P, N) chi² matrix.

        Returns (sums, maxes, upper_box_threshold, upper_threshold, pooled_null).
        """
        sums = np.sum(chi2_matrix, axis=1)
        maxes = np.max(chi2_matrix, axis=1)
        pooled_null = chi2_matrix.ravel().astype(np.float64)
        upper_box_threshold = float(np.percentile(pooled_null, 75))
        upper_threshold = float(np.percentile(pooled_null, 95))
        return sums, maxes, upper_box_threshold, upper_threshold, pooled_null

    def run(self, alignment_array, chi_square_calculator):
        """Run the permutation test and get chi-squared score percentiles and distribution."""
        print(f"Running {self.permutations} permutations.")
        print(f"Using {self.num_workers} worker(s) for permutation")
        chi2_matrix = self._permute_and_calculate_chi2(
            alignment_array, chi_square_calculator
        )
        return self._summarize_null(chi2_matrix)

    @property
    def _stratum_pools(self):
        """Per-taxon null pools from the last stratified run, or None."""
        if self._stratified_result is not None:
            return self._stratified_result.stratum_pools
        return None

    def compute_null(self, alignment_array, chi_square_calculator,
                     strategy="standard", alignment=None):
        """Compute the permutation null distribution.

        For standard: bulk column permutation (all taxa exchangeable).
        For similarity_stratified: within-stratum permutation with global
        baseline chi2.  Stores the StratifiedResult internally so that
        per-stratum pools are available via calc_empirical_pvalue().

        Parameters
        ----------
        alignment_array : (N, L) char array
        chi_square_calculator : ChiSquareCalculator
        strategy : str, "standard" or "similarity_stratified"
        alignment : Bio.Align.MultipleSeqAlignment, required for
            similarity_stratified (provides taxon names and raw sequences).

        Returns
        -------
        tuple : (sums, maxes, upper_box_threshold, upper_threshold, pooled_null)
        """
        if strategy == "similarity_stratified":
            from .stratified_permutation import run_similarity_stratified

            self._stratified_result = run_similarity_stratified(
                alignment_array,
                alignment,
                chi_square_calculator,
                permutation_test=self,
            )
            return self._stratified_result.as_standard_tuple()
        else:
            self._stratified_result = None
            return self.run(alignment_array, chi_square_calculator)

    def run_test(self, alignment_file, alignment_format, create_output=False,
                 strategy="standard"):
        """Run the permutation test on an alignment."""
        from .alignment_reader import AlignmentReader
        from .sequence_type_detector import SequenceTypeDetector
        from .chi_square_calculator import ChiSquareCalculator
        import time

        start_time = time.time()
        print(f"Permutation strategy: {strategy}")
        reader = AlignmentReader(alignment_file, alignment_format)
        alignment, alignment_array = reader.run()
        detector = SequenceTypeDetector()
        is_dna, char_set = detector.detect(alignment)
        self.alignment = alignment

        self.chi_square_calculator = ChiSquareCalculator(char_set, self.num_workers)
        sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2 = (
            self.compute_null(
                alignment_array, self.chi_square_calculator,
                strategy=strategy, alignment=alignment,
            )
        )

        from .utils import _robust_zscore

        row_counts = self.chi_square_calculator.calculate_row_counts(alignment_array)
        row_expected_observed = self.chi_square_calculator.calculate_expected_observed(
            row_counts
        )
        per_row_chi2 = self.chi_square_calculator.calculate_row_chi2(
            row_expected_observed, row_counts
        )

        # --- Stratification diagnostic (automatic for stratified) ---
        if (strategy == "similarity_stratified"
                and self._stratified_result is not None):
            from .stratification_diagnostic import (
                compare_null_distributions, print_diagnostic,
            )
            n_strata = self._stratified_result.diagnostics["n_strata_realizable"]
            strat_sums = self._stratified_result.sums
            observed_chi2 = float(np.sum(per_row_chi2))
            if n_strata > 1:
                print("Running standard permutation for diagnostic comparison.")
                pt_std = PermutationTest(self.num_workers, self.permutations)
                std_sums, _, _, _, _ = pt_std.run(
                    alignment_array, self.chi_square_calculator,
                )
                diag = compare_null_distributions(
                    observed_chi2, std_sums, strat_sums, n_strata,
                )
            else:
                p = float(np.sum(strat_sums >= observed_chi2) / len(strat_sums))
                med = float(np.median(strat_sums))
                diag = {
                    "valid": True, "warning": False,
                    "inflation_z": 0.0, "signal_consumed": 0.0,
                    "observed_chi2": observed_chi2,
                    "p_standard": p, "p_stratified": p,
                    "concordant": True,
                    "standard_median": med, "stratified_median": med,
                    "n_strata": 1, "alpha": 0.05,
                }
            print_diagnostic(diag)

        upper_chi_quantile = np.percentile(per_row_chi2, 95)

        median_z = _robust_zscore(np.mean(per_row_chi2), permutated_per_row_chi2)
        q95_z = _robust_zscore(upper_chi_quantile, permutated_per_row_chi2)

        print(
            f"Alignment chi2score: {(np.sum(per_row_chi2)):.2f} | "
            f"Permutations alignment chi2scores: {(min(sums)):.2f} - {(max(sums)):.2f} | "
            f"Empirical-P: {self.calc_empirical_pvalue(np.sum(per_row_chi2),sums)[0] }"
        )
        print(
            f"Alignment mean taxa chi2score: {(np.mean(per_row_chi2)):.2f} | "
            f"Permutations median taxa chi2score: {np.median(permutated_per_row_chi2):.2f} | "
            f"Mean Z-score: {median_z:.2f} | "
            f"q95 Z-score: {q95_z:.2f}"
        )
        # calculate zscores and empirical p-values for each row
        empirical_pvalues = self.calc_empirical_pvalue(
            per_row_chi2, permutated_per_row_chi2,
            per_taxon_pools=self._stratum_pools,
        )
        row_empirical_pvalue_dict = make_score_dict(
            per_row_chi2, permutated_per_row_chi2, empirical_pvalues, alignment,
            per_taxon_pools=self._stratum_pools,
        )
        # check for significant rows
        significant_list = [
            " ".join([str(t), str(row_empirical_pvalue_dict[t])])
            for t in row_empirical_pvalue_dict.keys()
            if row_empirical_pvalue_dict[t]["empirical_pvalue"] < 0.05
        ]

        print(
            f"Biased Taxa based on corrected empirical-P-values: {len(significant_list)} of {np.shape(per_row_chi2)[0]}"
        )
        if create_output:
            # remove file general file extension after last point (could be fasta or other)
            # form alignment file and add _score_dict.json
            output_tsv_file = alignment_file.replace(
                "." + alignment_file.split(".")[-1], "_scores.tsv"
            )
            write_score_dict_to_tsv(row_empirical_pvalue_dict, output_tsv_file)
            print(f"Printing taxa p-values and z-scores to file: {output_tsv_file}")
            # still print the top 5 z-score taxa to console
            print("Top 5 taxa z-scores and corresponding empirical p-values:")
            for t in list(row_empirical_pvalue_dict.keys())[:5]:
                print(
                    f"{t}\t{row_empirical_pvalue_dict[t]['zscore']}\t{row_empirical_pvalue_dict[t]['empirical_pvalue']}"
                )
        else:
            # print to console the zscore and pvalue per taxa from row_empirical_pvalue_dict
            print("Taxa z-scores and empirical p-values:")
            for t in row_empirical_pvalue_dict.keys():
                # print it in tabular format per taxon
                print(
                    f"{t}\t{row_empirical_pvalue_dict[t]['zscore']}\t{row_empirical_pvalue_dict[t]['empirical_pvalue']}"
                )
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Execution time for testing: {elapsed_time:.2f} seconds")
