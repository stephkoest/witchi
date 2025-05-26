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
            num_workers=1,
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
        self.tester = PermutationTest(num_workers=1, permutations=100)

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


if __name__ == "__main__":
    unittest.main()
