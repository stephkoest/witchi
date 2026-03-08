# tests/test_utils.py
import unittest
import os
import json
import tempfile
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from witchi.alignment_pruner.utils import (
    write_alignment,
    write_pruned_dict_to_tsv,
    write_score_dict_to_json,
    write_score_dict_to_tsv,
    _robust_zscore,
    make_score_dict,
)


class TestRobustZscore(unittest.TestCase):

    def test_typical_values(self):
        null = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        z = _robust_zscore(5.0, null)
        self.assertIsInstance(z, float)
        self.assertGreater(z, 0.0)

    def test_observed_equals_median_returns_zero(self):
        null = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        z = _robust_zscore(3.0, null)
        self.assertEqual(z, 0.0)

    def test_zero_mad_observed_equals_median(self):
        null = np.array([5.0, 5.0, 5.0, 5.0])
        z = _robust_zscore(5.0, null)
        self.assertEqual(z, 0.0)

    def test_zero_mad_observed_above_median(self):
        null = np.array([5.0, 5.0, 5.0, 5.0])
        z = _robust_zscore(6.0, null)
        self.assertEqual(z, np.inf)

    def test_zero_mad_observed_below_median(self):
        null = np.array([5.0, 5.0, 5.0, 5.0])
        z = _robust_zscore(4.0, null)
        self.assertEqual(z, -np.inf)

    def test_symmetric_null(self):
        null = np.arange(1, 102, dtype=float)  # odd length, median=51
        z_above = _robust_zscore(76.0, null)
        z_below = _robust_zscore(26.0, null)
        self.assertAlmostEqual(z_above, -z_below, places=10)


class TestMakeScoreDict(unittest.TestCase):

    def _make_alignment(self, names):
        records = [SeqRecord(Seq("ACGT"), id=n) for n in names]
        return MultipleSeqAlignment(records)

    def test_basic_output(self):
        names = ["tax1", "tax2", "tax3"]
        alignment = self._make_alignment(names)
        per_row_chi2 = np.array([10.0, 5.0, 1.0])
        null_pool = np.random.default_rng(0).uniform(0, 15, size=300)
        pvalues = [0.01, 0.5, 0.9]

        result = make_score_dict(per_row_chi2, null_pool, pvalues, alignment)
        self.assertEqual(set(result.keys()), set(names))
        for v in result.values():
            self.assertIn("empirical_pvalue", v)
            self.assertIn("zscore", v)
            self.assertNotIn("stratum", v)

    def test_sorted_by_zscore_descending(self):
        names = ["low", "mid", "high"]
        alignment = self._make_alignment(names)
        per_row_chi2 = np.array([1.0, 5.0, 20.0])
        null_pool = np.random.default_rng(0).uniform(0, 10, size=300)
        pvalues = [0.9, 0.5, 0.01]

        result = make_score_dict(per_row_chi2, null_pool, pvalues, alignment)
        zscores = [v["zscore"] for v in result.values()]
        self.assertEqual(zscores, sorted(zscores, reverse=True))

    def test_zscore_always_uses_pooled_null(self):
        names = ["a", "b"]
        alignment = self._make_alignment(names)
        per_row_chi2 = np.array([5.0, 10.0])
        pooled_null = np.random.default_rng(0).uniform(0, 10, size=200)
        pvalues = [0.1, 0.2]

        result = make_score_dict(
            per_row_chi2, pooled_null, pvalues, alignment,
        )
        self.assertEqual(len(result), 2)
        for v in result.values():
            self.assertIn("zscore", v)

    def test_with_stratum_info(self):
        names = ["x", "y"]
        alignment = self._make_alignment(names)
        per_row_chi2 = np.array([3.0, 7.0])
        null_pool = np.ones(100) * 5.0
        pvalues = [0.5, 0.5]
        name_to_stratum = {"x": 0, "y": 1}

        result = make_score_dict(
            per_row_chi2, null_pool, pvalues, alignment,
            name_to_stratum=name_to_stratum,
        )
        for name, v in result.items():
            self.assertIn("stratum", v)
            self.assertEqual(v["stratum"], name_to_stratum[name])


class TestWriteAlignment(unittest.TestCase):

    def test_write_and_read_back(self):
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1"),
            SeqRecord(Seq("ACGTACGT"), id="seq2"),
        ]
        alignment = MultipleSeqAlignment(records)
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            path = f.name
        try:
            write_alignment(alignment, path, format="fasta")
            self.assertTrue(os.path.exists(path))
            with open(path) as fh:
                content = fh.read()
            self.assertIn(">seq1", content)
            self.assertIn("ACGTACGT", content)
        finally:
            os.unlink(path)


class TestWritePrunedDictToTsv(unittest.TestCase):

    def test_writes_correct_headers_and_rows(self):
        prune_dict = {
            0: [1, 5, 100.0, 90.0, 10.0, 2],
            1: [2, 12, 95.0, 90.0, 5.0, 1],
        }
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False, mode="w") as f:
            path = f.name
        try:
            write_pruned_dict_to_tsv(prune_dict, path, "quartic")
            with open(path) as fh:
                lines = fh.readlines()
            self.assertIn("Iteration", lines[0])
            self.assertIn("Algorithm", lines[0])
            self.assertEqual(len(lines), 3)  # header + 2 rows
            self.assertIn("quartic", lines[1])
        finally:
            os.unlink(path)


class TestWriteScoreDictToJson(unittest.TestCase):

    def test_writes_valid_json_with_ndarrays(self):
        score_dict = {
            "before_real": np.array([1.0, 2.0, 3.0]),
            "after_real": np.array([0.5, 1.5, 2.5]),
        }
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            write_score_dict_to_json(score_dict, path)
            with open(path) as fh:
                loaded = json.load(fh)
            self.assertEqual(loaded["before_real"], [1.0, 2.0, 3.0])
        finally:
            os.unlink(path)


class TestWriteScoreDictToTsv(unittest.TestCase):

    def test_without_stratum(self):
        d = {
            "tax1": {"empirical_pvalue": 0.01, "zscore": 5.0},
            "tax2": {"empirical_pvalue": 0.5, "zscore": 1.0},
        }
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            path = f.name
        try:
            write_score_dict_to_tsv(d, path)
            with open(path) as fh:
                lines = fh.readlines()
            self.assertIn("Row", lines[0])
            self.assertNotIn("Stratum", lines[0])
            # sorted by descending abs zscore: tax1 (5.0) first
            self.assertTrue(lines[1].startswith("tax1"))
        finally:
            os.unlink(path)

    def test_with_stratum(self):
        d = {
            "tax1": {"empirical_pvalue": 0.01, "zscore": 3.0, "stratum": 0},
            "tax2": {"empirical_pvalue": 0.5, "zscore": 1.0, "stratum": 1},
        }
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            path = f.name
        try:
            write_score_dict_to_tsv(d, path)
            with open(path) as fh:
                lines = fh.readlines()
            self.assertIn("Stratum", lines[0])
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
