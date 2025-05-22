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
)


class AlignmentPruner:
    def __init__(
        self,
        file,
        format="fasta",
        max_residue_pruned=100,
        permutations=100,
        num_workers=2,
        top_n=10,
        pruning_algorithm="global",
        touchdown=False,
    ):
        self.file = file
        self.format = format
        self.max_residue_pruned = max_residue_pruned
        self.permutations = permutations
        self.num_workers = num_workers
        self.touchdown = touchdown
        self.top_n = top_n
        self.pruning_algorithm = pruning_algorithm
        self.chi_square_calculator = None

    def run(self):
        """Main method to run the recursive chi-square pruning."""
        start_time = time.time()
        reader = AlignmentReader(self.file, self.format)
        alignment, alignment_array = reader.run()

        is_dna, char_set = SequenceTypeDetector.detect(alignment)

        self.chi_square_calculator = ChiSquareCalculator(char_set, self.num_workers)
        self.permutation_test = PermutationTest(self.num_workers, self.permutations)
        sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2 = (
            self.permutation_test.run(alignment_array, self.chi_square_calculator)
        )
        median_perm_chi2 = np.median(permutated_per_row_chi2)
        mad_perm_chi2 = np.median(np.abs(permutated_per_row_chi2 - median_perm_chi2))
        pruned_alignment_array, prune_dict, score_dict = self.recursive_prune(
            alignment_array,
            sums,
            maxes,
            upper_box_threshold,
            upper_threshold,
            permutated_per_row_chi2,
            median_perm_chi2,
            mad_perm_chi2,
        )

        pruned_sequences = self.update_sequences(alignment, pruned_alignment_array)
        pruned_alignment = MultipleSeqAlignment(pruned_sequences)

        suffix = f"_{self.pruning_algorithm}_s{self.top_n}_pruned.fasta"
        output_alignment_pruned_file = os.path.splitext(self.file)[0] + suffix
        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", ".tsv"
        )
        write_alignment(pruned_alignment, output_alignment_pruned_file)
        write_pruned_dict_to_tsv(prune_dict, output_tsv_file, self.pruning_algorithm)

        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", "_score_dict.tsv"
        )
        # something to clean later, reintegrate into the main function
        empirical_pvalues = self.permutation_test.calc_empirical_pvalue(
            score_dict["after_real"], permutated_per_row_chi2
        )
        row_empirical_pvalue_dict = self.permutation_test.make_score_dict(
            score_dict["after_real"],
            permutated_per_row_chi2,
            empirical_pvalues,
            alignment,
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
        sums,
        permutated_per_row_chi2,
        upper_box_threshold,
        upper_threshold,
    ):
        if self.pruning_algorithm == "squared":
            initial_global_chi2 = self.chi_square_calculator.calculate_global_chi2(
                expected_observed, count_rows_array
            )
            chi2_differences = (
                self.chi_square_calculator.calculate_global_chi2_difference(
                    count_rows_array, alignment_array, initial_global_chi2
                )
            )
        elif self.pruning_algorithm == "wasserstein":
            wasserstein = self.chi_square_calculator.calculate_row_chi2_wasserstein(
                expected_observed, count_rows_array, permutated_per_row_chi2
            )
            chi2_differences = (
                self.chi_square_calculator.calculate_wasserstein_difference(
                    count_rows_array,
                    alignment_array,
                    wasserstein,
                    permutated_per_row_chi2,
                )
            )
            initial_global_chi2 = wasserstein
        elif self.pruning_algorithm == "quartic":
            initial_global_chi2 = (
                self.chi_square_calculator.calculate_quartic_row_global_chi2(
                    expected_observed, count_rows_array
                )
            )
            chi2_differences = (
                self.chi_square_calculator.calculate_quartic_chi2_difference(
                    count_rows_array, alignment_array, initial_global_chi2
                )
            )
        else:
            raise ValueError(f"Unknown pruning algorithm: {self.pruning_algorithm}")

        return initial_global_chi2, chi2_differences

    def recursive_prune(
        self,
        alignment_array,
        sums,
        maxes,
        upper_box_threshold,
        upper_threshold,
        permutated_per_row_chi2,
        median_perm_chi2,
        mad_perm_chi2,
    ):
        """Recursively prune the alignment array."""
        prune_dict = {}
        score_dict = {}
        iteration = 0
        removed_columns_count = 0
        alignment_size = alignment_array.shape[1]
        initial_topn = self.top_n
        original_indices = list(range(alignment_array.shape[1]))

        mean_perm_chi2 = np.mean(permutated_per_row_chi2)
        sd_perm_chi2 = np.std(permutated_per_row_chi2)

        while removed_columns_count < self.max_residue_pruned:
            stats = self._calculate_per_row_stats(alignment_array)
            count_rows_array = stats["count_rows"]
            expected_observed = stats["expected_observed"]
            per_row_chi2 = stats["per_row_chi2"]
            per_row_chi2_median = stats["median"]
            upper_chi_quantile = stats["q95"]
            global_chi2 = stats["sum"]

            if iteration == 0:
                score_dict["before_permuted"] = permutated_per_row_chi2
                score_dict["before_real"] = per_row_chi2

            empirical_pvalues = self.permutation_test.calc_empirical_pvalue(
                per_row_chi2, permutated_per_row_chi2
            )
            significant_count = sum(p <= 0.05 for p in empirical_pvalues)
            alignment_empiricalpvalue = self.permutation_test.calc_empirical_pvalue(
                np.sum(per_row_chi2), sums
            )[0]

            if per_row_chi2_median <= upper_box_threshold:
                if self.touchdown:
                    if self.top_n > 5:
                        self.top_n = 5

            if alignment_empiricalpvalue >= 0.95:
                print("Pruning complete. Exiting because of alignment p-value.")
                break

            initial_global_chi2, chi2_differences = self.prune(
                alignment_array,
                expected_observed,
                count_rows_array,
                sums,
                permutated_per_row_chi2,
                upper_box_threshold,
                upper_threshold,
            )
            # Sort the columns by the chi2 difference and get the top n columns to prune
            top_n_indices = np.argsort(list(chi2_differences.values()))[-self.top_n :]
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

            alignment_array = np.delete(alignment_array, top_n_indices, axis=1)

            original_indices = [
                i for j, i in enumerate(original_indices) if j not in top_n_indices
            ]

            print(
                f"Columns removed: {removed_columns_count}, {(removed_columns_count / alignment_size) * 100:.2f}% | "
                f"Biased taxa permutation: {significant_count} | "
                f"Mean z-score: {(np.mean(per_row_chi2) - mean_perm_chi2) / sd_perm_chi2:.2f} | "
                f"q95 z-score: {(upper_chi_quantile - upper_threshold) / (upper_threshold - mean_perm_chi2):.2f} | "
                f"Alignment p-value: {alignment_empiricalpvalue:.2f}"
            )

            if significant_count == 0:
                # quit when no significant taxa are left
                print("Pruning complete. Exiting because of taxa p-value.")
                self.top_n = initial_topn
                break
            iteration += 1

        score_dict["after_real"] = per_row_chi2

        return alignment_array, prune_dict, score_dict

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
