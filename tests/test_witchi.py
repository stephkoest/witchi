# tests/test_witchi.py
import unittest
from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.permutation_test import PermutationTest

class TestAlignmentPruner(unittest.TestCase):

    def setUp(self):
        self.pruner = AlignmentPruner(file='tests/data/example.nex', format='nexus', max_residue_pruned=100,
                                      permutations=100, num_workers=2, top_n=1, pruning_algorithm='wasserstein')

    def test_prune_alignment(self):
        pruned_alignment = self.pruner.run()

    def test_pruning_algorithm(self):
        self.pruner.pruning_algorithm = 'global'
        pruned_alignment = self.pruner.run()

class TestPermutationTest(unittest.TestCase):

    def setUp(self):
        self.tester = PermutationTest(num_workers=2, permutations=100)

    def test_run_test(self):
        result = self.tester.run_test(alignment_file='tests/data/example.nex', alignment_format='nexus')

    #def test_permutation_count(self): #TO FIX
    #    self.tester.permutations = 50
    #    result = self.tester.run_test(alignment_file='tests/data/example.nex', alignment_format='nexus')
    #    self.assertEqual(len(result['chi_square_scores']), 50)

if __name__ == '__main__':
    unittest.main()