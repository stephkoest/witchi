# tests/test_chi_square_calculator.py
import unittest
import numpy as np

from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator
from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
from witchi.alignment_pruner.alignment_reader import AlignmentReader


class TestChiSquareCalculator(unittest.TestCase):

    def setUp(self):
        reader = AlignmentReader("tests/data/example.nex", "nexus")
        self.alignment, self.alignment_array = reader.run()
        _, self.char_set = SequenceTypeDetector.detect(self.alignment)
        self.calc = ChiSquareCalculator(self.char_set, num_workers=1)

    def test_calculate_row_counts_shape(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        n_chars = len(self.char_set)
        n_taxa = self.alignment_array.shape[0]
        self.assertEqual(counts.shape, (n_chars, n_taxa))

    def test_calculate_row_counts_positive(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        self.assertTrue(np.all(counts > 0))

    def test_calculate_expected_observed_shape(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        self.assertEqual(expected.shape, counts.shape)

    def test_calculate_expected_observed_positive(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        self.assertTrue(np.all(expected > 0))

    def test_global_chi2_scalar_positive(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        chi2 = self.calc.calculate_global_chi2(expected, counts)
        self.assertIsInstance(float(chi2), float)
        self.assertGreater(chi2, 0.0)

    def test_row_chi2_shape_and_nonnegative(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        row_chi2 = self.calc.calculate_row_chi2(expected, counts)
        n_taxa = self.alignment_array.shape[0]
        self.assertEqual(row_chi2.shape, (n_taxa,))
        self.assertTrue(np.all(row_chi2 >= 0))

    def test_row_chi2_sums_to_global(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        global_chi2 = self.calc.calculate_global_chi2(expected, counts)
        row_chi2 = self.calc.calculate_row_chi2(expected, counts)
        np.testing.assert_allclose(np.sum(row_chi2), global_chi2, rtol=1e-10)

    def test_quartic_row_global_chi2_is_sum_of_squares(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        row_chi2 = self.calc.calculate_row_chi2(expected, counts)
        quartic = self.calc.calculate_quartic_row_global_chi2(expected, counts)
        np.testing.assert_allclose(quartic, np.sum(row_chi2 ** 2), rtol=1e-10)

    def test_wasserstein_nonnegative(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        row_chi2 = self.calc.calculate_row_chi2(expected, counts)
        K = 50
        null_quantiles = np.quantile(row_chi2, np.linspace(0, 1, K + 2)[1:-1])
        w = self.calc.calculate_row_chi2_wasserstein(expected, counts, null_quantiles)
        self.assertGreaterEqual(w, 0.0)

    def test_global_chi2_difference_returns_dict(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        global_chi2 = self.calc.calculate_global_chi2(expected, counts)
        diffs = self.calc.calculate_global_chi2_difference(
            counts, self.alignment_array, global_chi2
        )
        self.assertIsInstance(diffs, dict)
        self.assertEqual(len(diffs), self.alignment_array.shape[1])

    def test_quartic_chi2_difference_returns_dict(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        quartic = self.calc.calculate_quartic_row_global_chi2(expected, counts)
        diffs = self.calc.calculate_quartic_chi2_difference(
            counts, self.alignment_array, quartic
        )
        self.assertIsInstance(diffs, dict)
        self.assertEqual(len(diffs), self.alignment_array.shape[1])

    def test_wasserstein_difference_returns_dict(self):
        counts = self.calc.calculate_row_counts(self.alignment_array)
        expected = self.calc.calculate_expected_observed(counts)
        row_chi2 = self.calc.calculate_row_chi2(expected, counts)
        K = 50
        null_quantiles = np.quantile(row_chi2, np.linspace(0, 1, K + 2)[1:-1])
        w = self.calc.calculate_row_chi2_wasserstein(expected, counts, null_quantiles)
        diffs = self.calc.calculate_wasserstein_difference(
            counts, self.alignment_array, w, null_quantiles
        )
        self.assertIsInstance(diffs, dict)
        self.assertEqual(len(diffs), self.alignment_array.shape[1])


class TestChiSquareCalculatorSynthetic(unittest.TestCase):
    """Tests with synthetic data for edge cases."""

    def test_uniform_composition_low_chi2(self):
        char_set = np.array(list("ACGT"))
        calc = ChiSquareCalculator(char_set, num_workers=1)
        rng = np.random.default_rng(42)
        # 4 taxa, 100 columns, roughly uniform composition
        alignment_array = rng.choice(list("ACGT"), size=(4, 100))
        counts = calc.calculate_row_counts(alignment_array)
        expected = calc.calculate_expected_observed(counts)
        row_chi2 = calc.calculate_row_chi2(expected, counts)
        # uniform composition should yield low chi2 per taxon
        self.assertTrue(np.all(row_chi2 < 50))

    def test_biased_taxon_has_highest_chi2(self):
        char_set = np.array(list("ACGT"))
        calc = ChiSquareCalculator(char_set, num_workers=1)
        rng = np.random.default_rng(42)
        # 3 uniform taxa + 1 GC-biased taxon
        uniform = rng.choice(list("ACGT"), size=(3, 200))
        biased = rng.choice(list("GC"), size=(1, 200))
        alignment_array = np.vstack([uniform, biased])
        counts = calc.calculate_row_counts(alignment_array)
        expected = calc.calculate_expected_observed(counts)
        row_chi2 = calc.calculate_row_chi2(expected, counts)
        self.assertEqual(np.argmax(row_chi2), 3)


if __name__ == "__main__":
    unittest.main()
