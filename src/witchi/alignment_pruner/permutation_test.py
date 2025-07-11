import numpy as np
from joblib import Parallel, delayed

from .utils import write_score_dict_to_tsv, make_score_dict


class PermutationTest:
    def __init__(self, num_workers_permute, permutations):
        self.is_dna = None
        self.chi_square_calculator = None
        self.num_workers = num_workers_permute
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
        print(f"Using {self.num_workers} worker(s) for permutation")
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
            f"Mean Z-score: {(np.mean(per_row_chi2) - mean_perm_chi2) / sd_perm_chi2:.2f} | "
            f"q95 Z-score: {(upper_chi_quantile - upper_threshold) / sd_perm_chi2:.2f}"
        )
        # calculate zscores and empirical p-values for each row
        empirical_pvalues = self.calc_empirical_pvalue(
            per_row_chi2, permutated_per_row_chi2
        )
        row_empirical_pvalue_dict = make_score_dict(
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
        # write_score_dict_to_json(sorted_row_chi2, "row_chi2_scores.json")
