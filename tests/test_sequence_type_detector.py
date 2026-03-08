# tests/test_sequence_type_detector.py
import unittest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from witchi.alignment_pruner.sequence_type_detector import SequenceTypeDetector


class TestSequenceTypeDetector(unittest.TestCase):

    def _make_alignment(self, seqs):
        records = [SeqRecord(Seq(s), id=f"seq{i}") for i, s in enumerate(seqs)]
        return MultipleSeqAlignment(records)

    def test_detect_dna(self):
        alignment = self._make_alignment(["ACGTACGT", "GCTAGCTA"])
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertTrue(is_dna)
        np.testing.assert_array_equal(sorted(char_set), sorted(SequenceTypeDetector.DNA_BASES))

    def test_detect_protein(self):
        alignment = self._make_alignment(["MVLSPADKTNVK", "MWLSGEDKSNIK"])
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertFalse(is_dna)
        np.testing.assert_array_equal(
            sorted(char_set), sorted(SequenceTypeDetector.AMINO_ACIDS)
        )

    def test_dna_with_gaps_and_ambiguity(self):
        alignment = self._make_alignment(["ACGT-NNN", "GC?TACGT"])
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertTrue(is_dna)

    def test_lowercase_detected_as_dna(self):
        alignment = self._make_alignment(["acgtacgt", "gctagcta"])
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertTrue(is_dna)

    def test_all_ambiguous_raises(self):
        alignment = self._make_alignment(["NNN---", "???..."])
        with self.assertRaises(ValueError):
            SequenceTypeDetector.detect(alignment)

    def test_unexpected_characters_default_to_protein(self):
        # 'J' is not in DNA or standard amino acids
        alignment = self._make_alignment(["ACGTJKLM", "MVLSJKLM"])
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertFalse(is_dna)

    def test_from_real_dna_file(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader

        reader = AlignmentReader("tests/data/example.nex", "nexus")
        alignment, _ = reader.run()
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertTrue(is_dna)

    def test_from_real_protein_file(self):
        from witchi.alignment_pruner.alignment_reader import AlignmentReader

        reader = AlignmentReader("tests/data/example_protein.fasta", "fasta")
        alignment, _ = reader.run()
        is_dna, char_set = SequenceTypeDetector.detect(alignment)
        self.assertFalse(is_dna)


if __name__ == "__main__":
    unittest.main()
