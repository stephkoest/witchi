# tests/test_witchi.py
import unittest
from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.permutation_test import PermutationTest
import os
import numpy as np


class TestAlignmentPruner(unittest.TestCase):

    def setUp(self):
        self.pruner = AlignmentPruner(
            file="tests/data/example.nex",
            format="nexus",
            max_residue_pruned=100,
            permutations=100,
            num_workers_chisq=1,
            num_workers_permute=1,
            top_n=1,
            pruning_algorithm="quartic",
        )

    def test_prune_alignment(self):
        self.pruner.run()
        base = "tests/data/example_quartic_s1_pruned"
        assert os.path.exists(base + ".fasta")
        assert os.path.exists(base + ".tsv")
        assert os.path.exists(base + "_scores.tsv")
        assert os.path.exists(base + "_score_dict.json")

    def test_pruning_algorithm(self):
        self.pruner.pruning_algorithm = "squared"
        self.pruner.run()

    def test_recursive_prune_returns_expected_keys(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader
        from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
        from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator

        reader = AlignmentReader(self.pruner.file, self.pruner.format)
        alignment, alignment_array = reader.run()
        _, char_set = SequenceTypeDetector.detect(alignment)
        self.pruner.chi_square_calculator = ChiSquareCalculator(char_set, 1)
        self.pruner.permutation_test = PermutationTest(1, 10)

        sums, _, _, _, perm_chi2 = self.pruner.permutation_test.run(
            alignment_array, self.pruner.chi_square_calculator
        )

        result = self.pruner.recursive_prune(alignment_array, sums, 99.0, perm_chi2)
        assert isinstance(result, tuple)
        assert len(result) == 3
        assert "after_real" in result[2]

    def test_invalid_algorithm_raises(self):
        self.pruner.pruning_algorithm = "invalid"
        with self.assertRaises(ValueError):
            self.pruner.run()


class TestPermutationTest(unittest.TestCase):

    def setUp(self):
        self.tester = PermutationTest(num_workers_permute=1, permutations=100)

    def test_run_test(self):
        self.tester.run_test(
            alignment_file="tests/data/example.nex", alignment_format="nexus"
        )

    def test_run_output_structure(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader
        from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
        from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator

        reader = AlignmentReader("tests/data/example.nex", "nexus")
        alignment, alignment_array = reader.run()
        _, char_set = SequenceTypeDetector.detect(alignment)
        calc = ChiSquareCalculator(char_set, num_workers=1)

        sums, maxes, q75, q95, flat_chi2 = self.tester.run(alignment_array, calc)

        self.assertEqual(sums.ndim, 1)
        self.assertEqual(maxes.ndim, 1)
        self.assertIsInstance(q75, float)
        self.assertIsInstance(q95, float)
        self.assertEqual(flat_chi2.ndim, 1)

    def test_calc_empirical_pvalue_vector_and_scalar(self):
        observed = np.array([5, 10, 15])
        null = np.array([1, 2, 5, 10, 20, 30])
        pvals = self.tester.calc_empirical_pvalue(observed, null)
        self.assertEqual(len(pvals), 3)
        self.assertTrue(all(0 <= p for p in pvals))

        single = self.tester.calc_empirical_pvalue(12, null)
        self.assertTrue(0 <= single[0] <= 1)

    def test_run_with_few_permutations(self):
        self.tester.permutations = 2
        from witchi.alignment_pruner.alignment_reader import AlignmentReader
        from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
        from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator

        alignment, array = AlignmentReader("tests/data/example.nex", "nexus").run()
        _, char_set = SequenceTypeDetector.detect(alignment)
        calc = ChiSquareCalculator(char_set, 1)
        output = self.tester.run(array, calc)
        self.assertEqual(len(output), 5)


class TestStratifiedPermutation(unittest.TestCase):
    """Tests for the similarity-stratified permutation pipeline."""

    def setUp(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader
        from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
        from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator

        reader = AlignmentReader("tests/data/example.nex", "nexus")
        self.alignment, self.alignment_array = reader.run()
        _, self.char_set = SequenceTypeDetector.detect(self.alignment)
        self.chi_calc = ChiSquareCalculator(self.char_set, num_workers=1)
        self.N = self.alignment_array.shape[0]  # 5 taxa
        self.P = 50  # keep tests fast

    # --- msa_treecut_stratification ---

    def test_msa_strata_returns_all_taxa(self):
        from witchi.alignment_pruner.msa_treecut_stratification import msa_strata

        names = [rec.id for rec in self.alignment]
        raw_seqs = [str(rec.seq) for rec in self.alignment]
        result = msa_strata(raw_seqs, names, min_stratum_size=2)
        self.assertEqual(set(result.keys()), set(names))
        # stratum indices are non-negative integers
        self.assertTrue(all(v >= 0 for v in result.values()))

    def test_msa_strata_diagnostics(self):
        from witchi.alignment_pruner.msa_treecut_stratification import msa_strata

        names = [rec.id for rec in self.alignment]
        raw_seqs = [str(rec.seq) for rec in self.alignment]
        result, diag = msa_strata(raw_seqs, names, min_stratum_size=2,
                                  return_diagnostics=True)
        self.assertIn("n_strata", diag)
        self.assertIn("n_strata_natural", diag)
        self.assertIn("isolation", diag)
        self.assertEqual(len(diag["isolation"]), self.N)

    def test_msa_strata_respects_min_stratum_size(self):
        from witchi.alignment_pruner.msa_treecut_stratification import msa_strata

        names = [rec.id for rec in self.alignment]
        raw_seqs = [str(rec.seq) for rec in self.alignment]
        # With min_stratum_size >= N, everything must land in one stratum
        result = msa_strata(raw_seqs, names, min_stratum_size=self.N)
        self.assertEqual(len(set(result.values())), 1)

    # --- stratified_permutation core ---

    def _make_stratified_result(self):
        """Helper: run stratified permutation and return StratifiedResult."""
        from witchi.alignment_pruner.stratified_permutation import (
            run_similarity_stratified,
        )
        pt = PermutationTest(num_workers_permute=1, permutations=self.P)
        return run_similarity_stratified(
            self.alignment_array, self.alignment, self.chi_calc,
            permutation_test=pt,
        )

    def test_run_similarity_stratified_returns_stratified_result(self):
        from witchi.alignment_pruner.stratified_permutation import StratifiedResult
        result = self._make_stratified_result()
        self.assertIsInstance(result, StratifiedResult)

    def test_stratified_result_standard_tuple_shape(self):
        result = self._make_stratified_result()
        sums, maxes, q75, q95, pooled = result.as_standard_tuple()

        self.assertEqual(sums.shape, (self.P,))
        self.assertEqual(maxes.shape, (self.P,))
        self.assertIsInstance(q75, float)
        self.assertIsInstance(q95, float)
        self.assertEqual(pooled.ndim, 1)
        self.assertEqual(len(pooled), self.P * self.N)

    def test_stratified_result_chi2_matrix_shape(self):
        result = self._make_stratified_result()
        self.assertEqual(result.chi2_matrix.shape, (self.P, self.N))

    def test_stratified_result_stratum_pools_cover_all_taxa(self):
        result = self._make_stratified_result()
        self.assertEqual(set(result.stratum_pools.keys()),
                         set(range(self.N)))
        for pool in result.stratum_pools.values():
            self.assertGreater(len(pool), 0)

    def test_stratified_result_bin_ids_shape(self):
        result = self._make_stratified_result()
        self.assertEqual(len(result.bin_ids), self.N)
        self.assertTrue(np.all(result.bin_ids >= 0))

    # --- per-stratum p-value ---

    def test_calc_empirical_pvalue_with_per_taxon_pools(self):
        result = self._make_stratified_result()
        # compute observed chi2
        rows = self.chi_calc.calculate_row_counts(self.alignment_array)
        exp = self.chi_calc.calculate_expected_observed(rows)
        obs = self.chi_calc.calculate_row_chi2(exp, rows)

        tester = PermutationTest(num_workers_permute=1, permutations=self.P)
        pooled = result.as_standard_tuple()[4]
        pvals = tester.calc_empirical_pvalue(
            obs, pooled, per_taxon_pools=result.stratum_pools,
        )
        self.assertEqual(len(pvals), self.N)
        for p in pvals:
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    # --- integration through PermutationTest dispatch ---

    def test_compute_null_stratified_returns_5tuple(self):
        tester = PermutationTest(num_workers_permute=1, permutations=self.P)
        result = tester.compute_null(
            self.alignment_array, self.chi_calc,
            strategy="similarity_stratified", alignment=self.alignment,
        )
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 5)

    def test_compute_null_stores_stratified_result(self):
        tester = PermutationTest(num_workers_permute=1, permutations=self.P)
        tester.compute_null(
            self.alignment_array, self.chi_calc,
            strategy="similarity_stratified", alignment=self.alignment,
        )
        self.assertIsNotNone(tester._stratified_result)

    def test_calc_empirical_pvalue_uses_strata_after_compute_null(self):
        tester = PermutationTest(num_workers_permute=1, permutations=self.P)
        _, _, _, _, pooled = tester.compute_null(
            self.alignment_array, self.chi_calc,
            strategy="similarity_stratified", alignment=self.alignment,
        )
        rows = self.chi_calc.calculate_row_counts(self.alignment_array)
        exp = self.chi_calc.calculate_expected_observed(rows)
        obs = self.chi_calc.calculate_row_chi2(exp, rows)

        pvals = tester.calc_empirical_pvalue(
            obs, pooled, per_taxon_pools=tester._stratum_pools,
        )
        self.assertEqual(len(pvals), self.N)
        for p in pvals:
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    def test_standard_vs_stratified_both_produce_valid_pvalues(self):
        """Both strategies should produce valid p-values on the same data."""
        tester_std = PermutationTest(num_workers_permute=1, permutations=self.P)
        tester_strat = PermutationTest(num_workers_permute=1, permutations=self.P)

        _, _, _, _, pooled_std = tester_std.compute_null(
            self.alignment_array, self.chi_calc, strategy="standard",
        )
        _, _, _, _, pooled_strat = tester_strat.compute_null(
            self.alignment_array, self.chi_calc,
            strategy="similarity_stratified", alignment=self.alignment,
        )

        rows = self.chi_calc.calculate_row_counts(self.alignment_array)
        exp = self.chi_calc.calculate_expected_observed(rows)
        obs = self.chi_calc.calculate_row_chi2(exp, rows)

        pvals_std = tester_std.calc_empirical_pvalue(
            obs, pooled_std, per_taxon_pools=tester_std._stratum_pools,
        )
        pvals_strat = tester_strat.calc_empirical_pvalue(
            obs, pooled_strat, per_taxon_pools=tester_strat._stratum_pools,
        )

        for pvals in (pvals_std, pvals_strat):
            self.assertEqual(len(pvals), self.N)
            for p in pvals:
                self.assertGreaterEqual(p, 0.0)
                self.assertLessEqual(p, 1.0)

    # --- full CLI-level integration ---

    def test_run_test_stratified(self):
        """run_test with similarity_stratified completes without error."""
        tester = PermutationTest(num_workers_permute=1, permutations=self.P)
        tester.run_test(
            alignment_file="tests/data/example.nex",
            alignment_format="nexus",
            strategy="similarity_stratified",
        )

    def test_prune_stratified(self):
        """Full pruning with similarity_stratified completes and writes outputs."""
        pruner = AlignmentPruner(
            file="tests/data/example.nex",
            format="nexus",
            max_residue_pruned=10,
            permutations=self.P,
            num_workers_chisq=1,
            num_workers_permute=1,
            top_n=1,
            pruning_algorithm="quartic",
            strategy="similarity_stratified",
        )
        pruner.run()
        base = "tests/data/example_quartic_s1_similarity_stratified_pruned"
        self.assertTrue(os.path.exists(base + ".fasta"))
        self.assertTrue(os.path.exists(base + ".tsv"))


class TestFreedmanLaneChi2(unittest.TestCase):
    """Unit tests for the _freedman_lane_chi2 helper."""

    def test_single_stratum_equals_standard_chi2(self):
        """With one stratum containing all taxa, FL should equal standard chi2."""
        from witchi.alignment_pruner.permutation_test import _freedman_lane_chi2

        rng = np.random.default_rng(99)
        count_rows_array = rng.integers(5, 50, size=(4, 6)).astype(np.float64)

        strata = [np.array([0, 1, 2, 3, 4, 5])]
        fl_chi2 = _freedman_lane_chi2(count_rows_array, strata)

        # Standard chi2
        total_per_char = count_rows_array.sum(axis=1)
        global_freq = total_per_char / total_per_char.sum()
        taxon_totals = count_rows_array.sum(axis=0)
        expected = global_freq[:, np.newaxis] * taxon_totals[np.newaxis, :]
        std_chi2 = ((count_rows_array - expected) ** 2 / expected).sum(axis=0)

        np.testing.assert_allclose(fl_chi2, std_chi2, rtol=1e-10)

    def test_two_strata_reduces_chi2(self):
        """Compositionally distinct strata should yield FL chi2 <= standard."""
        from witchi.alignment_pruner.permutation_test import _freedman_lane_chi2

        # Stratum 0: GC-rich, Stratum 1: AT-rich
        count_rows_array = np.array(
            [
                [10, 12, 40, 38, 42],
                [40, 42, 10, 12, 11],
                [38, 40, 12, 10, 13],
                [12, 10, 38, 40, 34],
            ],
            dtype=np.float64,
        )

        strata = [np.array([0, 1]), np.array([2, 3, 4])]
        fl_chi2 = _freedman_lane_chi2(count_rows_array, strata)

        # Standard chi2
        total_per_char = count_rows_array.sum(axis=1)
        global_freq = total_per_char / total_per_char.sum()
        taxon_totals = count_rows_array.sum(axis=0)
        expected = global_freq[:, np.newaxis] * taxon_totals[np.newaxis, :]
        std_chi2 = ((count_rows_array - expected) ** 2 / expected).sum(axis=0)

        self.assertTrue(np.all(fl_chi2 <= std_chi2 + 1e-10))
        self.assertTrue(np.any(fl_chi2 < std_chi2 - 1e-6))

    def test_output_shape_and_nonnegative(self):
        """Output should be (N,) and all non-negative."""
        from witchi.alignment_pruner.permutation_test import _freedman_lane_chi2

        rng = np.random.default_rng(7)
        count_rows_array = rng.integers(1, 30, size=(4, 10)).astype(np.float64)
        strata = [np.array([0, 1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])]
        fl_chi2 = _freedman_lane_chi2(count_rows_array, strata)

        self.assertEqual(fl_chi2.shape, (10,))
        self.assertTrue(np.all(fl_chi2 >= 0))


class TestStratificationDiagnostic(unittest.TestCase):
    """Tests for the stratification validity diagnostic."""

    def setUp(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader
        from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector
        from witchi.alignment_pruner.chi_square_calculator import ChiSquareCalculator

        reader = AlignmentReader("tests/data/example.nex", "nexus")
        self.alignment, self.alignment_array = reader.run()
        _, self.char_set = SequenceTypeDetector.detect(self.alignment)
        self.chi_calc = ChiSquareCalculator(self.char_set, num_workers=1)
        self.N = self.alignment_array.shape[0]

    def test_diagnostic_returns_expected_keys(self):
        """Return dict has all documented keys with correct types."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=20, num_workers=1,
        )
        expected_keys = {
            "valid", "warning", "inflation_z", "signal_consumed",
            "observed_chi2", "p_standard", "p_stratified", "concordant",
            "standard_median", "stratified_median",
            "n_strata", "trivial", "alpha",
        }
        self.assertEqual(set(result.keys()), expected_keys)
        self.assertIsInstance(result["valid"], bool)
        self.assertIsInstance(result["warning"], bool)
        self.assertIsInstance(result["concordant"], bool)
        self.assertIsInstance(result["trivial"], bool)
        self.assertIsInstance(result["observed_chi2"], float)
        self.assertIsInstance(result["inflation_z"], float)
        self.assertIsInstance(result["signal_consumed"], float)

    def test_diagnostic_p_values_in_range(self):
        """Both p-values are in [0, 1]."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=20, num_workers=1,
        )
        for key in ("p_standard", "p_stratified"):
            self.assertGreaterEqual(result[key], 0.0)
            self.assertLessEqual(result[key], 1.0)

    def test_diagnostic_trivial_with_few_taxa(self):
        """With 5 taxa, likely 1 stratum -> trivial=True, valid=True."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=20, num_workers=1,
        )
        if result["n_strata"] <= 1:
            self.assertTrue(result["trivial"])
            self.assertTrue(result["valid"])
            self.assertFalse(result["warning"])
            self.assertTrue(result["concordant"])
            self.assertEqual(result["inflation_z"], 0.0)
            self.assertEqual(result["p_standard"], result["p_stratified"])

    def test_diagnostic_concordance_logic(self):
        """valid == concordant, and concordance matches p-value agreement."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=50, num_workers=1,
        )
        alpha = result["alpha"]
        expected_concordant = (
            (result["p_standard"] < alpha) == (result["p_stratified"] < alpha)
        )
        self.assertEqual(result["concordant"], expected_concordant)
        self.assertEqual(result["valid"], result["concordant"])

    def test_diagnostic_observed_chi2_is_positive(self):
        """Observed chi2 should be a positive number."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=20, num_workers=1,
        )
        self.assertGreater(result["observed_chi2"], 0.0)

    def test_diagnostic_warning_requires_both_criteria(self):
        """warning requires inflation_z >= 3 AND signal_consumed >= 0.10."""
        from witchi.alignment_pruner.stratification_diagnostic import (
            diagnose_stratification_validity,
        )
        result = diagnose_stratification_validity(
            self.alignment_array, self.alignment, self.chi_calc,
            permutations=50, num_workers=1,
        )
        expected_warning = (
            result["inflation_z"] >= 3.0 and result["signal_consumed"] >= 0.10
        )
        self.assertEqual(result["warning"], expected_warning)
        self.assertGreaterEqual(result["inflation_z"], 0.0)
        self.assertGreaterEqual(result["signal_consumed"], 0.0)


if __name__ == "__main__":
    unittest.main()
