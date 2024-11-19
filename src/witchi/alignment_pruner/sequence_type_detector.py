import numpy as np

class SequenceTypeDetector:
    DNA_BASES = np.array(list("ACGT"))
    AMINO_ACIDS = np.array(list("ACDEFGHIKLMNPQRSTVWY"))

    def __init__(self):
        self.is_dna = None
        self.char_set = None

    def detect_sequence_type(self, alignment):
        """Detect whether the alignment is DNA or Protein based on sequence characters."""
        chars_to_remove = {'N', 'n', '-', '?', 'X', 'x', '.'}
        seq_set = set()

        for record in alignment:
            seq_set.update(list(str(record.seq.upper())))

        for char in chars_to_remove:
            seq_set.discard(char)

        if seq_set.issubset(self.DNA_BASES):
            self.is_dna = True
            self.char_set = self.DNA_BASES
            print("Detected DNA sequence.")
        else:
            self.is_dna = False
            self.char_set = self.AMINO_ACIDS
            print("Detected Protein sequence.")