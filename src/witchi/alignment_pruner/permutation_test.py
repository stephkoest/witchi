import json
import numpy as np
from joblib import Parallel, delayed
import csv


class PermutationTest:
    def __init__(self, num_workers, permutations):
        self.is_dna = None
        self.chi_square_calculator = None
        self.num_workers = num_workers
        self.permutations = permutations

    def calc_empirical_pvalue(self, per_row_chi2, permutated_per_row_chi2):
        """Calculate empirical p-value."""
        empirical_p_list = []
        # check if array, indicating per taxon chi2 scores
        if isinstance(per_row_chi2, np.ndarray):
            for i in range(len(per_row_chi2)):
                # divide number of permutated_per_row_chi2 larger than per_row_chi2[i]
                # by the number of taxa in permutations to get the probability of the chi2 score
                empirical_p = (
                    np.sum(permutated_per_row_chi2 >= per_row_chi2[i])
                    / len(permutated_per_row_chi2)
                ) * len(per_row_chi2)
                empirical_p_list.append(empirical_p)
        # if not array, indicating total alignment chi2 scores
        else:
            empirical_p = np.sum(permutated_per_row_chi2 >= per_row_chi2) / len(
                permutated_per_row_chi2
            )
            empirical_p_list.append(empirical_p)

        return empirical_p_list

    def _permute_and_calculate_chi2(self, alignment_array, chi_square_calculator):
        def permute_single_calculate_chi2(i):
            iter_seed = 12345 + i
            rng = np.random.default_rng(iter_seed)
            permuted_array = np.apply_along_axis(rng.permutation, 0, alignment_array)
            count_rows_array = chi_square_calculator.calculate_row_counts(
                permuted_array
            )
            expected_observed = chi_square_calculator.calculate_expected_observed(
                count_rows_array
            )
            per_row_chi2 = chi_square_calculator.calculate_row_chi2(
                expected_observed, count_rows_array
            )
            return per_row_chi2

        results = Parallel(n_jobs=self.num_workers)(
            delayed(permute_single_calculate_chi2)(i) for i in range(self.permutations)
        )
        return np.array(results)

    def run(self, alignment_array, chi_square_calculator):
        """Run the permutation test and get chi-squared score percentiles and distribution."""
        print(f"Running {self.permutations} permutations.")
        permutated_per_row_chi2 = self._permute_and_calculate_chi2(
            alignment_array, chi_square_calculator
        )

        maxes = np.max(permutated_per_row_chi2, axis=1)
        # gett sums for every permutation
        sums = np.sum(permutated_per_row_chi2, axis=1)
        # Flatten the list of chi-squared scores
        permutated_per_row_chi2 = np.concatenate(permutated_per_row_chi2)

        upper_box_threshold = np.percentile(permutated_per_row_chi2, 75)
        upper_threshold = np.percentile(permutated_per_row_chi2, 95)

        return (
            sums,
            maxes,
            upper_box_threshold,
            upper_threshold,
            permutated_per_row_chi2,
        )

    def run_test(self, alignment_file, alignment_format, create_output=False):
        """Run the permutation test on an alignment."""
        from .alignment_reader import AlignmentReader
        from .sequence_type_detector import SequenceTypeDetector
        from .chi_square_calculator import ChiSquareCalculator
        import time

        start_time = time.time()
        reader = AlignmentReader(alignment_file, alignment_format)
        alignment, alignment_array = reader.run()
        detector = SequenceTypeDetector()
        is_dna, char_set = detector.detect(alignment)
        self.alignment = alignment

        self.chi_square_calculator = ChiSquareCalculator(char_set, self.num_workers)
        sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2 = (
            self.run(alignment_array, self.chi_square_calculator)
        )
        mean_perm_chi2 = np.mean(permutated_per_row_chi2)
        sd_perm_chi2 = np.std(permutated_per_row_chi2)

        row_counts = self.chi_square_calculator.calculate_row_counts(alignment_array)
        row_expected_observed = self.chi_square_calculator.calculate_expected_observed(
            row_counts
        )
        per_row_chi2 = self.chi_square_calculator.calculate_row_chi2(
            row_expected_observed, row_counts
        )

        upper_chi_quantile = np.percentile(per_row_chi2, 95)

        print(
            f"Alignment chi2score: {(np.sum(per_row_chi2)):.2f} | "
            f"Permutations alignment chi2scores: {(min(sums)):.2f} - {(max(sums)):.2f} | "
            f"Empirical-P: {self.calc_empirical_pvalue(np.sum(per_row_chi2),sums)[0] }"
        )
        print(
            f"Alignment mean taxa chi2score: {(np.mean(per_row_chi2)):.2f} | "
            f"Permutations mean taxa chi2score: {(mean_perm_chi2):.2f} | "
            f"Mean z-score: {(np.mean(per_row_chi2) - mean_perm_chi2) / sd_perm_chi2:.2f} | "
            f"q95 z-score: {(upper_chi_quantile - upper_threshold) / (upper_threshold - mean_perm_chi2):.2f}"
        )
        # calculate zscores and empirical p-values for each row
        empirical_pvalues = self.calc_empirical_pvalue(
            per_row_chi2, permutated_per_row_chi2
        )
        row_empirical_pvalue_dict = self.make_score_dict(
            per_row_chi2, permutated_per_row_chi2, empirical_pvalues, alignment
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
            self.write_score_dict_to_tsv(row_empirical_pvalue_dict, output_tsv_file)
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
        # self.write_score_dict_to_json(sorted_row_chi2, "row_chi2_scores.json")

    def make_score_dict(
        self, per_row_chi2, permutated_per_row_chi2, empirical_pvalues, alignment
    ):
        """Make a dictionary of chi-squared scores for each row in the alignment."""
        # Extract row names
        per_row_chi2 = np.array(per_row_chi2)
        row_names = [record.id for record in alignment]
        mean_perm_chi2 = np.mean(permutated_per_row_chi2)
        sd_perm_chi2 = np.std(permutated_per_row_chi2)
        # Sort the dictionary by chi-squared scores in descending order
        zscores = (per_row_chi2 - mean_perm_chi2) / sd_perm_chi2

        # Calculate pvalues for each row
        row_empirical_pvalue_dict = {
            row_names[i]: {
                "empirical_pvalue": empirical_pvalues[i],
                "zscore": zscores[i],
            }
            for i in range(len(row_names))
        }
        # sort row_empirical_pvalue_dict by z-score in decending order
        row_empirical_pvalue_dict = dict(
            sorted(
                row_empirical_pvalue_dict.items(),
                key=lambda item: item[1]["zscore"],
                reverse=True,
            )
        )

        return row_empirical_pvalue_dict

    def write_score_dict_to_tsv(self, dictionary, file_name):
        """Write the score dictionary to a TSV file, ordered by descending absolute z-score."""
        # Sort the dictionary by descending absolute z-score
        sorted_dict = dict(
            sorted(
                dictionary.items(),
                key=lambda item: abs(item[1]["zscore"]),
                reverse=True,
            )
        )

        with open(file_name, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t")
            writer.writerow(["Row", "Empirical-Pvalue", "Z-Score"])
            for row, values in sorted_dict.items():
                writer.writerow([row, values["empirical_pvalue"], values["zscore"]])

    def write_score_dict_to_json(self, dictionary, file_name):
        """Write the score dictionary to a JSON file."""

        def convert_ndarray(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

        with open(file_name, "w") as jsonfile:
            json.dump(dictionary, jsonfile, indent=4, default=convert_ndarray)
