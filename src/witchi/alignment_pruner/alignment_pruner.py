import numpy as np
import os
import time
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

from .chi_square_calculator import ChiSquareCalculator
from .permutation_test import PermutationTest
from .sequence_type_detector import SequenceTypeDetector
from .alignment_reader import AlignmentReader
from .utils import (
    write_alignment,
    write_pruned_dict_to_tsv,
    write_score_dict_to_json,
    write_score_dict_to_tsv,
    make_score_dict,
)


class AlignmentPruner:
    _DELTA_NULL_PERMUTATIONS = 20

    def __init__(
        self,
        file,
        format="fasta",
        max_residue_pruned=100,
        permutations=100,
        num_workers_chisq=1,
        num_workers_permute=1,
        top_n=10,
        pruning_algorithm="squared",
        touchdown=False,
        strategy="standard",
        strict=False,
        delta_null=True,
    ):
        self.file = file
        self.format = format
        self.max_residue_pruned = max_residue_pruned
        self.permutations = permutations
        self.num_workers_chisq = num_workers_chisq
        self.num_workers_permute = num_workers_permute
        self.touchdown = touchdown
        self.strict = strict
        self.top_n = top_n
        self.pruning_algorithm = pruning_algorithm
        self.strategy = strategy
        self.delta_null = delta_null
        self.chi_square_calculator = None
        self.alignment_size = None
        self.initial_top_n = self.top_n
        self._null_max_deltas = None
        self.prune_strategies = {
            "squared": self._prune_squared,
            "quartic": self._prune_quartic,
            "wasserstein": self._prune_wasserstein,
        }

    def run(self):
        """Main method to run the recursive chi-square pruning."""
        start_time = time.time()
        print("Starting pruning workflow")
        print(f"Chi² calculator using {self.num_workers_chisq} worker(s)")
        print(f"Permutation test using {self.num_workers_permute} worker(s)")
        print(f"Pruning algorithm: {self.pruning_algorithm}")
        print(f"Permutation strategy: {self.strategy}")
        print(f"Initial top_n: {self.top_n}")
        reader = AlignmentReader(self.file, self.format)
        alignment, alignment_array = reader.run()

        is_dna, char_set = SequenceTypeDetector.detect(alignment)

        self.chi_square_calculator = ChiSquareCalculator(
            char_set, self.num_workers_chisq
        )
        self.permutation_test = PermutationTest(
            self.num_workers_permute, self.permutations
        )
        (
            sums,
            maxes,
            upper_box_threshold,
            upper_threshold,
            permutated_per_row_chi2,
        ) = self.permutation_test.compute_null(
            alignment_array,
            self.chi_square_calculator,
            strategy=self.strategy,
            alignment=alignment,
        )

        # Precompute Z-score parameters from pooled null (consistent with _robust_zscore)
        self._null_z_mean = float(np.mean(permutated_per_row_chi2))
        _median = float(np.median(permutated_per_row_chi2))
        _mad = float(np.median(np.abs(permutated_per_row_chi2 - _median)))
        self._null_z_scale = _mad / 0.6745

        # Compress null Z-scores to quantiles for Wasserstein pruning.
        # For stratified strategy, target the stratified null (tree-aware shape);
        # for standard, target the standard null.
        if (
            self.strategy == "similarity_stratified"
            and self.permutation_test._stratified_result is not None
        ):
            target_null = self.permutation_test._stratified_result.pooled_null
        else:
            target_null = permutated_per_row_chi2
        null_z = (target_null - self._null_z_mean) / self._null_z_scale
        K = min(200, len(null_z))
        positions = np.linspace(0, 1, K + 2)[1:-1]
        self._null_z_quantiles = np.quantile(null_z, positions)

        if self.delta_null:
            print("Computing null delta distribution for stopping criterion...")
            self._compute_null_deltas(alignment_array)

        pruned_alignment_array, prune_dict, score_dict = self.recursive_prune(
            alignment_array,
            sums,
            upper_threshold,
            permutated_per_row_chi2,
            alignment=alignment,
        )

        pruned_sequences = self.update_sequences(alignment, pruned_alignment_array)
        pruned_alignment = MultipleSeqAlignment(pruned_sequences)
        # build output filename suffix
        parts = [f"_{self.pruning_algorithm}_s{self.top_n}"]
        if self.strategy != "standard":
            parts.append(f"_{self.strategy}")
        if self.touchdown:
            parts.append("_touchdown")
        parts.append("_pruned.fasta")
        suffix = "".join(parts)
        output_alignment_pruned_file = os.path.splitext(self.file)[0] + suffix
        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", ".tsv"
        )
        write_alignment(pruned_alignment, output_alignment_pruned_file)
        write_pruned_dict_to_tsv(prune_dict, output_tsv_file, self.pruning_algorithm)

        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", ".tsv"
        )
        empirical_pvalues = self.permutation_test.calc_empirical_pvalue(
            score_dict["after_real"],
            permutated_per_row_chi2,
            per_taxon_pools=self.permutation_test._stratum_pools,
        )
        name_to_stratum = None
        if self.permutation_test._stratified_result is not None:
            name_to_stratum = self.permutation_test._stratified_result.diagnostics.get(
                "name_to_stratum"
            )
        row_empirical_pvalue_dict = make_score_dict(
            score_dict["after_real"],
            permutated_per_row_chi2,
            empirical_pvalues,
            alignment,
            name_to_stratum=name_to_stratum,
        )
        output_score_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", "_scores.tsv"
        )
        write_score_dict_to_tsv(row_empirical_pvalue_dict, output_score_tsv_file)
        # NOW MAKE JSON
        output_json_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", "_score_dict.json"
        )
        write_score_dict_to_json(score_dict, output_json_file)
        print(f"Pruned {len(prune_dict.keys())} columns.")
        print(f"Pruned alignment saved to {output_alignment_pruned_file}")
        print(f"Columns pruned in order saved to: {output_tsv_file}")
        print(f"Taxa p-values and z-scores printed to {output_score_tsv_file}")
        print(f"Chiscore values from permutation testing printed to {output_json_file}")

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.2f} seconds")

    def prune(
        self,
        alignment_array,
        expected_observed,
        count_rows_array,
        permutated_per_row_chi2,
    ):
        try:
            return self.prune_strategies[self.pruning_algorithm](
                alignment_array,
                expected_observed,
                count_rows_array,
                permutated_per_row_chi2,
            )
        except KeyError:
            raise ValueError(f"Unknown pruning algorithm: {self.pruning_algorithm}")

    def _prune_squared(self, alignment_array, expected_observed, count_rows_array, *_):
        initial_global_chi2 = self.chi_square_calculator.calculate_global_chi2(
            expected_observed, count_rows_array
        )
        chi2_differences = self.chi_square_calculator.calculate_global_chi2_difference(
            count_rows_array, alignment_array, initial_global_chi2
        )
        return initial_global_chi2, chi2_differences

    def _prune_quartic(self, alignment_array, expected_observed, count_rows_array, *_):
        initial_global_chi2 = (
            self.chi_square_calculator.calculate_quartic_row_global_chi2(
                expected_observed, count_rows_array
            )
        )
        chi2_differences = self.chi_square_calculator.calculate_quartic_chi2_difference(
            count_rows_array, alignment_array, initial_global_chi2
        )
        return initial_global_chi2, chi2_differences

    def _prune_wasserstein(
        self,
        alignment_array,
        expected_observed,
        count_rows_array,
        permutated_per_row_chi2,
    ):
        wasserstein = self.chi_square_calculator.calculate_row_zscore_wasserstein(
            expected_observed,
            count_rows_array,
            self._null_z_quantiles,
            self._null_z_mean,
            self._null_z_scale,
        )
        chi2_differences = (
            self.chi_square_calculator.calculate_wasserstein_zscore_difference(
                count_rows_array,
                alignment_array,
                wasserstein,
                self._null_z_quantiles,
                self._null_z_mean,
                self._null_z_scale,
            )
        )
        return wasserstein, chi2_differences

    def _compute_null_deltas(self, alignment_array):
        """Compute null distribution of max per-column delta scores.

        For each of P_null permuted alignments, compute per-column deltas
        using the active pruning algorithm, then record the maximum delta.
        The resulting distribution tests whether observed deltas during
        pruning are distinguishable from chance.
        """
        p_null = self._DELTA_NULL_PERMUTATIONS
        max_deltas = np.empty(p_null, dtype=np.float64)

        # Determine strata for stratified permutation
        strata_indices = None
        sr = self.permutation_test._stratified_result
        if self.strategy == "similarity_stratified" and sr is not None:
            bin_ids = sr.bin_ids
            n_strata = int(np.max(bin_ids)) + 1
            strata_indices = [np.where(bin_ids == k)[0] for k in range(n_strata)]
            strata_indices = [s for s in strata_indices if len(s) >= 2]

        for i in range(p_null):
            # Fresh seeds, offset past the existing permutation seeds
            iter_seed = 12345 + self.permutations + i
            rng = np.random.default_rng(iter_seed)

            if strata_indices is None:
                permuted_array = np.apply_along_axis(
                    rng.permutation, 0, alignment_array
                )
            else:
                permuted_array = alignment_array.copy()
                for s_indices in strata_indices:
                    sub = permuted_array[s_indices, :]
                    permuted_array[s_indices, :] = np.apply_along_axis(
                        rng.permutation, 0, sub
                    )

            count_rows = self.chi_square_calculator.calculate_row_counts(permuted_array)
            expected = self.chi_square_calculator.calculate_expected_observed(
                count_rows
            )

            if self.pruning_algorithm == "squared":
                initial = self.chi_square_calculator.calculate_global_chi2(
                    expected, count_rows
                )
                deltas = self.chi_square_calculator.calculate_global_chi2_difference(
                    count_rows, permuted_array, initial
                )
            elif self.pruning_algorithm == "quartic":
                initial = self.chi_square_calculator.calculate_quartic_row_global_chi2(
                    expected, count_rows
                )
                deltas = self.chi_square_calculator.calculate_quartic_chi2_difference(
                    count_rows, permuted_array, initial
                )
            elif self.pruning_algorithm == "wasserstein":
                initial = self.chi_square_calculator.calculate_row_zscore_wasserstein(
                    expected,
                    count_rows,
                    self._null_z_quantiles,
                    self._null_z_mean,
                    self._null_z_scale,
                )
                deltas = (
                    self.chi_square_calculator.calculate_wasserstein_zscore_difference(
                        count_rows,
                        permuted_array,
                        initial,
                        self._null_z_quantiles,
                        self._null_z_mean,
                        self._null_z_scale,
                    )
                )
            else:
                raise ValueError(f"Unknown pruning algorithm: {self.pruning_algorithm}")

            max_deltas[i] = max(deltas.values())
            print(
                f"  Null delta permutation {i + 1}/{p_null}: "
                f"max_delta={max_deltas[i]:.6f}"
            )

        self._null_max_deltas = max_deltas
        print(
            f"Null max-delta distribution (n={p_null}): "
            f"mean={np.mean(max_deltas):.6f}, "
            f"p95={np.percentile(max_deltas, 95):.6f}"
        )

    def recursive_prune(
        self,
        alignment_array,
        sums,
        upper_threshold,
        permutated_per_row_chi2,
        alignment=None,
    ):
        """Recursively prune the alignment array."""
        prune_dict = {}
        score_dict = {}
        iteration = 0
        removed_columns_count = 0
        self.alignment_size = alignment_array.shape[1]
        original_indices = list(range(alignment_array.shape[1]))

        from .utils import _robust_zscore

        top_n_indices = []
        initial_global_chi2 = None
        chi2_differences = {}

        while removed_columns_count < self.max_residue_pruned:
            stats = self._calculate_per_row_stats(alignment_array)
            count_rows_array = stats["count_rows"]
            expected_observed = stats["expected_observed"]
            per_row_chi2 = stats["per_row_chi2"]
            upper_chi_quantile = stats["q95"]
            global_chi2 = stats["sum"]

            alignment_empirical_p, significant_count = self._calc_empirical_pvals(
                stats, sums, permutated_per_row_chi2
            )

            if iteration == 0:
                score_dict["before_permuted"] = permutated_per_row_chi2
                score_dict["before_real"] = per_row_chi2
                sr = self.permutation_test._stratified_result
                if sr is not None:
                    score_dict["stratified_permuted"] = sr.pooled_null
                    score_dict["stratification"] = {
                        "taxon_names": [rec.id for rec in alignment],
                        "bin_ids": sr.bin_ids.tolist(),
                        "n_strata": int(sr.diagnostics["n_strata_realizable"]),
                    }
            else:
                # only write to prune_dict after the first pruning iteration
                for col in top_n_indices:
                    pruned_col = original_indices[col]
                    prune_dict[pruned_col] = [
                        iteration,
                        pruned_col,
                        global_chi2,
                        initial_global_chi2,
                        chi2_differences[col],
                        significant_count,
                    ]
                    removed_columns_count += 1

                original_indices = [
                    i for j, i in enumerate(original_indices) if j not in top_n_indices
                ]

            should_stop, stop_reason = self._check_stopping_criteria(
                stats, permutated_per_row_chi2, alignment_empirical_p, significant_count
            )

            mean_z = _robust_zscore(np.mean(per_row_chi2), permutated_per_row_chi2)
            q95_z = _robust_zscore(upper_chi_quantile, permutated_per_row_chi2)

            print(
                f"Columns removed: {removed_columns_count}, "
                f"{(removed_columns_count / self.alignment_size) * 100:.2f}% | "
                f"Biased taxa permutation: {significant_count} | "
                f"Mean Z-score: {mean_z:.2f} | "
                f"q95 Z-score: {q95_z:.2f} | "
                f"Alignment p-value: {alignment_empirical_p:.2f}"
            )

            if should_stop:
                print(f"Pruning complete. Exiting because of {stop_reason}.")
                break

            initial_global_chi2, chi2_differences = self.prune(
                alignment_array,
                expected_observed,
                count_rows_array,
                permutated_per_row_chi2,
            )

            # Sort candidates in descending delta order, take up to top_n
            sorted_candidates = sorted(
                chi2_differences.items(), key=lambda x: x[1], reverse=True
            )[: self.top_n]

            # Delta null: walk candidates and remove those above the null max.
            # Stops at first failure within the batch, giving graceful
            # touchdown — batch size naturally shrinks as we approach
            # convergence. Exits the pruning loop only when rank-1 itself
            # fails (zero columns justified in this iteration).
            if self._null_max_deltas is not None:
                justified = []
                last_p = None
                for col, delta in sorted_candidates:
                    last_p = np.sum(self._null_max_deltas >= delta) / len(
                        self._null_max_deltas
                    )
                    if last_p > 0.05:
                        break
                    justified.append(col)

                if not justified:
                    max_delta = sorted_candidates[0][1]
                    print(
                        f"Pruning complete. Exiting because of delta not "
                        f"significant (p={last_p:.3f}, "
                        f"max_delta={max_delta:.6f})."
                    )
                    break

                if len(justified) < len(sorted_candidates):
                    fail_delta = sorted_candidates[len(justified)][1]
                    print(
                        f"  Partial batch: removing {len(justified)}/"
                        f"{len(sorted_candidates)} columns "
                        f"(rank {len(justified) + 1} fails at "
                        f"p={last_p:.3f}, delta={fail_delta:.6f})"
                    )
            else:
                justified = [col for col, _ in sorted_candidates]

            top_n_indices = np.array(justified)
            alignment_array = np.delete(alignment_array, top_n_indices, axis=1)
            iteration += 1

        score_dict["after_real"] = per_row_chi2
        self.top_n = self.initial_top_n
        return alignment_array, prune_dict, score_dict

    def _check_stopping_criteria(
        self, stats, permuted_chi2, alignment_empirical_p, significant_count
    ):
        """
        Check whether pruning should stop based on empirical p-values.
        Returns: (should_stop: bool, reason: str, significant_count: int)
        """
        if self.touchdown:
            if stats["median"] <= np.percentile(permuted_chi2, 99):
                new_top_n = int(self.alignment_size / 1000)
                if new_top_n < self.top_n:
                    self.top_n = new_top_n

        if alignment_empirical_p > 0.05 and significant_count == 0:
            return True, "convergence"

        if not self.strict:
            if alignment_empirical_p > 0.05:
                return True, "alignment p-value"

        if significant_count == 0:
            return True, "no significant taxa"

        return False, None

    def _calc_empirical_pvals(self, stats, permuted_sums, permuted_chi2):
        """Calculate empirical p-values for alignment and taxa."""
        alignment_empirical_p = self.permutation_test.calc_empirical_pvalue(
            stats["sum"], permuted_sums
        )[0]
        empirical_pvals = self.permutation_test.calc_empirical_pvalue(
            stats["per_row_chi2"],
            permuted_chi2,
            per_taxon_pools=self.permutation_test._stratum_pools,
        )
        significant_count = sum(p <= 0.05 for p in empirical_pvals)
        return alignment_empirical_p, significant_count

    def _calculate_per_row_stats(self, alignment_array):
        """Calculate row-wise chi² statistics and summary metrics."""
        count_rows_array = self.chi_square_calculator.calculate_row_counts(
            alignment_array
        )
        expected_observed = self.chi_square_calculator.calculate_expected_observed(
            count_rows_array
        )
        per_row_chi2 = self.chi_square_calculator.calculate_row_chi2(
            expected_observed, count_rows_array
        )

        return {
            "count_rows": count_rows_array,
            "expected_observed": expected_observed,
            "per_row_chi2": per_row_chi2,
            "median": np.median(per_row_chi2),
            "mean": np.mean(per_row_chi2),
            "q95": np.percentile(per_row_chi2, 95),
            "sum": np.sum(per_row_chi2),
        }

    def update_sequences(self, alignment, pruned_alignment_array):
        """Update sequences with pruned alignments."""
        pruned_sequences = []

        for i, record in enumerate(alignment):
            pruned_seq = "".join(pruned_alignment_array[i, :])
            record.seq = Seq(pruned_seq)
            pruned_sequences.append(record)

        return pruned_sequences
