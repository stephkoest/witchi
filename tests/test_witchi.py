# tests/test_witchi.py
import unittest
from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.permutation_test import PermutationTest


class TestAlignmentPruner(unittest.TestCase):

    def setUp(self):
        self.pruner = AlignmentPruner(
            file="tests/data/example.nex",
            format="nexus",
            max_residue_pruned=100,
            permutations=100,
            num_workers=2,
            top_n=1,
            pruning_algorithm="quartic",
        )

    def test_prune_alignment(self):
        self.pruner.run()

    def test_pruning_algorithm(self):
        self.pruner.pruning_algorithm = "squared"
        self.pruner.run()


class TestPermutationTest(unittest.TestCase):

    def setUp(self):
        self.tester = PermutationTest(num_workers=2, permutations=100)

    def test_run_test(self):
        self.tester.run_test(
            alignment_file="tests/data/example.nex", alignment_format="nexus"
        )


if __name__ == "__main__":
    unittest.main()
