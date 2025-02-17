import json
import numpy as np
from joblib import Parallel, delayed
from Bio.Align import MultipleSeqAlignment

class PermutationTest:
    def __init__(self, num_workers, permutations):
        self.is_dna = None
        self.chi_square_calculator = None
        self.num_workers = num_workers
        self.permutations = permutations

    def calc_pseudo_pvalue(self, per_row_chi2, permutated_per_row_chi2):
        """Calculate pseudo p-value."""
        pseudo_p_list = []
        #check if array, indicating per taxon chi2 scores
        if isinstance(per_row_chi2, np.ndarray):
            for i in range(len(per_row_chi2)):
                #divide number of permutated_per_row_chi2 larger than per_row_chi2[i] by the number of taxa in permutations to get the probability of the chi2 score
                pseudo_p = (np.sum(permutated_per_row_chi2 >= per_row_chi2[i]) / len(permutated_per_row_chi2)) * len(per_row_chi2)
                pseudo_p_list.append(pseudo_p)
        #if not array, indicating total alignment chi2 scores
        else:
            pseudo_p = (np.sum(permutated_per_row_chi2 >= per_row_chi2) / len(permutated_per_row_chi2))
            pseudo_p_list.append(pseudo_p)

        return pseudo_p_list

    def run(self, alignment_array, chi_square_calculator):
        """Run the permutation test and get chi-squared score percentiles and distribution."""
        def permute_single_calculate_chi2(i):
            iter_seed = 12345 + i
            rng = np.random.default_rng(iter_seed)
            permuted_array = np.apply_along_axis(rng.permutation, 0, alignment_array)
            count_rows_array = chi_square_calculator.calculate_row_counts(permuted_array)
            expected_observed = chi_square_calculator.calculate_expected_observed(count_rows_array)
            per_row_chi2 = chi_square_calculator.calculate_row_chi2(expected_observed, count_rows_array)
            return per_row_chi2

        print(f"Running {self.permutations} permutations.")
        permutated_per_row_chi2 = Parallel(n_jobs=self.num_workers)(
            delayed(permute_single_calculate_chi2)(i) for i in range(self.permutations)
        )

        permutated_per_row_chi2 = np.array(permutated_per_row_chi2)
        maxes = np.max(permutated_per_row_chi2, axis=1)
        #gett sums for every permutation
        sums = np.sum(permutated_per_row_chi2, axis=1)
        # Flatten the list of chi-squared scores
        permutated_per_row_chi2 = np.concatenate(permutated_per_row_chi2)

        upper_box_threshold = np.percentile(permutated_per_row_chi2, 75)
        upper_threshold = np.percentile(permutated_per_row_chi2, 95)

        return sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2

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
        detector.detect_sequence_type(alignment)
        self.alignment = alignment
        self.is_dna = detector.is_dna
        char_set = detector.char_set
        self.chi_square_calculator = ChiSquareCalculator(char_set, self.num_workers)
        sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2 = self.run(alignment_array, self.chi_square_calculator)
        mean_perm_chi2 = np.mean(permutated_per_row_chi2)
        sd_perm_chi2 = np.std(permutated_per_row_chi2)

        row_counts = self.chi_square_calculator.calculate_row_counts(alignment_array)
        row_expected_observed = self.chi_square_calculator.calculate_expected_observed(row_counts)
        per_row_chi2 = self.chi_square_calculator.calculate_row_chi2(row_expected_observed, row_counts)
        significant_count_permutation = np.sum(per_row_chi2 > upper_threshold)

        upper_chi_quantile = np.percentile(per_row_chi2, 95)

        # Extract row names
        row_names = [record.id for record in alignment]
        # Sort the dictionary by chi-squared scores in descending order
        #sorted_row_chi2 = dict(sorted(row_chi2_dict.items(), key=lambda item: item[1], reverse=True))
        pseudo_pvalues = self.calc_pseudo_pvalue(per_row_chi2, permutated_per_row_chi2)
        zscores = (per_row_chi2 - mean_perm_chi2) / sd_perm_chi2
        print(f"Mean z-score: {(np.mean(per_row_chi2) - mean_perm_chi2) / sd_perm_chi2:.2f} | q95 z-score: {(upper_chi_quantile - upper_threshold) / (upper_threshold - mean_perm_chi2):.2f}")
        print(f"Permutations alignment chi2scores: {(min(sums)):.2f} - {(max(sums)):.2f} | Alignment chi2score: {(np.sum(per_row_chi2)):.2f} | Pseudo-P: {self.calc_pseudo_pvalue(np.sum(per_row_chi2),sums)[0] }")
        print(f"Permutations mean chi2score: {(mean_perm_chi2):.2f} | Alignment mean chi2score: {(np.mean(per_row_chi2)):.2f} ")
        # Calculate pvalues for each row
        row_pseudo_pvalue_dict = {row_names[i]: {'pseudo_pvalue': pseudo_pvalues[i], 'zscore': zscores[i]} for i in
                                  range(len(row_names))}
        #check for significant rows
        significant_list = [" ".join([str(t), str(row_pseudo_pvalue_dict[t])]) for t in row_pseudo_pvalue_dict.keys() if row_pseudo_pvalue_dict[t]['pseudo_pvalue'] < 0.05]
        significant_string = ", ".join(significant_list)
        print(f"Biased Taxa based on corrected pseudo-p-values: {len(significant_list)} of {np.shape(per_row_chi2)[0]}")
        if create_output:
            #remove file general file extension after last point (could be fasta or other) form alignment file and add _score_dict.json
            output_json_file = alignment_file.replace('.' + alignment_file.split('.')[-1], '_score_dict.json')
            self.write_score_dict_to_json(row_pseudo_pvalue_dict, output_json_file)
            print(f"Printing taxa p-values and z-scores to file: {output_json_file}")
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Execution time for testing: {elapsed_time:.2f} sec