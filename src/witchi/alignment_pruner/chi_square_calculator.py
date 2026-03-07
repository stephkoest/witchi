import numpy as np
from joblib import Parallel, delayed


class ChiSquareCalculator:
    def __init__(self, char_set, num_workers=2):
        self.char_set = char_set
        self.num_workers = num_workers

    def calculate_row_counts(self, alignment_array):
        """
        Count occurrences of all characters in the provided char_set for each column
        of the alignment array, ignoring gaps ('-').
        """
        count_rows_array = np.array(
            [(alignment_array == char).sum(axis=1) for char in self.char_set]
        )
        if np.any(count_rows_array == 0):
            # Avoid division by zero
            count_rows_array = count_rows_array + 1

        return count_rows_array

    def calculate_expected_observed(self, count_rows_array):
        """
        Fully vectorized calculation of expected and observed frequencies for each column in
        the alignment using precomputed character counts.
        """
        total_chars = np.sum(count_rows_array, axis=1)
        total_chars_sum = np.sum(total_chars)
        global_frequencies = total_chars / total_chars_sum
        column_totals = np.sum(count_rows_array, axis=0, keepdims=True)
        expected_frequencies = global_frequencies[:, np.newaxis] * column_totals
        if np.any(expected_frequencies == 0):
            # Avoid division by zero
            expected_frequencies = expected_frequencies + 1

        return expected_frequencies

    def calculate_global_chi2(self, expected_values, count_rows_array):
        """Calculate the global chi-squared score for the entire alignment."""
        # Calculate chi-squared values
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        total_chi2 = np.sum(chi2_values)

        return total_chi2

    def calculate_row_chi2(self, expected_values, count_rows_array):
        """Calculate the chi-squared score for all taxa."""
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        chi2_sum = np.sum(chi2_values, axis=0)
        return chi2_sum

    def calculate_row_chi2_outlyingness_sum(
        self, expected_values, count_rows_array, median_perm_chi2, mad_perm_chi2
    ):
        """Calculate the chi-squared score for all taxa."""
        per_row_chi2 = self.calculate_row_chi2(expected_values, count_rows_array)
        outlyingness_sum = np.sum(
            self.calculate_perm_outlyingness(
                per_row_chi2, median_perm_chi2, mad_perm_chi2
            )
        )

        return outlyingness_sum

    def calculate_row_chi2_wasserstein(
        self, expected_values, count_rows_array, null_quantiles
    ):
        """Wasserstein-1 distance between observed per-taxon chi² and null.

        Parameters
        ----------
        null_quantiles : (K,) array
            Pre-compressed quantiles of the null distribution at K
            uniformly-spaced positions.
        """
        per_row_chi2 = self.calculate_row_chi2(expected_values, count_rows_array)
        K = len(null_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        obs_quantiles = np.quantile(per_row_chi2, positions)
        return float(np.mean(np.abs(obs_quantiles - null_quantiles)))

    def calculate_quartic_row_global_chi2(self, expected_values, count_rows_array):
        """Calculate the chi-squared score for all taxa."""
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        chi2_sum = np.sum(chi2_values, axis=0)
        quartic_chi2_sum = np.sum(chi2_sum**2)

        return quartic_chi2_sum

    def calculate_global_chi2_difference(
        self, count_rows_array, alignment_array, initial_global_chi2
    ):
        """Calculate the chi-squared difference for each column based on its impact on the overall chi-squared score."""

        def compute_global_chi2_difference(col_idx):
            col_count = count_rows_array - self.calculate_row_counts(
                alignment_array[:, col_idx, np.newaxis]
            )
            return self.calculate_global_chi2(
                self.calculate_expected_observed(col_count), col_count
            )

        chi2_differences = initial_global_chi2 - np.array(
            Parallel(n_jobs=self.num_workers)(
                delayed(compute_global_chi2_difference)(col_idx)
                for col_idx in range(alignment_array.shape[1])
            )
        )

        return dict(enumerate(chi2_differences))

    def calculate_outlyingness_difference(
        self,
        count_rows_array,
        alignment_array,
        outlyingness_sum,
        median_perm_chi2,
        mad_perm_chi2,
    ):
        """
        Calculate the chi-squared difference for each column based on its impact
        on the chi-squared score of significant rows.
        """

        def compute_row_chi2_outlyingness_difference(col_idx):
            col_count = count_rows_array - self.calculate_row_counts(
                alignment_array[:, col_idx, np.newaxis]
            )
            return self.calculate_row_chi2_outlyingness_sum(
                self.calculate_expected_observed(col_count),
                col_count,
                median_perm_chi2,
                mad_perm_chi2,
            )

        chi2_differences = outlyingness_sum - np.array(
            Parallel(n_jobs=self.num_workers)(
                delayed(compute_row_chi2_outlyingness_difference)(col_idx)
                for col_idx in range(alignment_array.shape[1])
            )
        )

        return dict(enumerate(chi2_differences))

    def calculate_wasserstein_difference(
        self, count_rows_array, alignment_array, wasserstein, null_quantiles
    ):
        """Calculate the wasserstein difference for each column based on its impact on the chi-squared score of rows."""

        def compute_row_chi2_wasserstein_difference(col_idx):
            col_count = count_rows_array - self.calculate_row_counts(
                alignment_array[:, col_idx, np.newaxis]
            )
            return self.calculate_row_chi2_wasserstein(
                self.calculate_expected_observed(col_count),
                col_count,
                null_quantiles,
            )

        chi2_differences = wasserstein - np.array(
            Parallel(n_jobs=self.num_workers)(
                delayed(compute_row_chi2_wasserstein_difference)(col_idx)
                for col_idx in range(alignment_array.shape[1])
            )
        )

        return dict(enumerate(chi2_differences))

    def calculate_row_zscore_wasserstein(
        self, expected_values, count_rows_array, null_z_quantiles, null_z_mean, null_z_scale
    ):
        """Wasserstein-1 distance between observed per-taxon Z-scores and null Z distribution.

        Converts per-taxon chi² to Z-scores using the pooled null parameters,
        then compares quantiles against the null Z quantiles.
        """
        per_row_chi2 = self.calculate_row_chi2(expected_values, count_rows_array)
        z_scores = (per_row_chi2 - null_z_mean) / null_z_scale
        K = len(null_z_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        obs_z_quantiles = np.quantile(z_scores, positions)
        return float(np.mean(np.abs(obs_z_quantiles - null_z_quantiles)))

    def calculate_wasserstein_zscore_difference(
        self, count_rows_array, alignment_array, wasserstein, null_z_quantiles, null_z_mean, null_z_scale
    ):
        """Calculate the wasserstein Z-score difference for each column."""

        def compute_row_zscore_wasserstein_difference(col_idx):
            col_count = count_rows_array - self.calculate_row_counts(
                alignment_array[:, col_idx, np.newaxis]
            )
            return self.calculate_row_zscore_wasserstein(
                self.calculate_expected_observed(col_count),
                col_count,
                null_z_quantiles,
                null_z_mean,
                null_z_scale,
            )

        chi2_differences = wasserstein - np.array(
            Parallel(n_jobs=self.num_workers)(
                delayed(compute_row_zscore_wasserstein_difference)(col_idx)
                for col_idx in range(alignment_array.shape[1])
            )
        )

        return dict(enumerate(chi2_differences))

    def calculate_quartic_chi2_difference(
        self, count_rows_array, alignment_array, initial_global_chi2
    ):
        """
        Calculate the chi-squared difference for each column
        based on its impact on the squared chi-squared score of rows.
        """

        def compute_quartic_difference(col_idx):
            col_count = count_rows_array - self.calculate_row_counts(
                alignment_array[:, col_idx, np.newaxis]
            )
            return self.calculate_quartic_row_global_chi2(
                self.calculate_expected_observed(col_count), col_count
            )

        chi2_differences = initial_global_chi2 - np.array(
            Parallel(n_jobs=self.num_workers)(
                delayed(compute_quartic_difference)(col_idx)
                for col_idx in range(alignment_array.shape[1])
            )
        )

        return dict(enumerate(chi2_differences))
