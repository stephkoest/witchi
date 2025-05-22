import numpy as np


class SequenceTypeDetector:
    DNA_BASES = np.array(list("ACGT"))
    AMINO_ACIDS = np.array(list("ACDEFGHIKLMNPQRSTVWY"))

    @staticmethod
    def detect(alignment):
        chars_to_remove = {"N", "n", "-", "?", "X", "x", "."}
        seq_set = set()

        for record in alignment:
            seq_set.update(list(str(record.seq.upper())))

        seq_set -= chars_to_remove

        if not seq_set:
            print(
                "⚠️  Could not detect sequence type: all residues were ambiguous or removed."
            )
            raise ValueError("No valid sequence characters detected.")

        if seq_set.issubset(SequenceTypeDetector.DNA_BASES):
            print("Detected DNA sequence.")
            return True, SequenceTypeDetector.DNA_BASES

        elif seq_set.issubset(SequenceTypeDetector.AMINO_ACIDS):
            print("Detected amino acid sequence.")
            return False, SequenceTypeDetector.AMINO_ACIDS

        else:
            print(f"⚠️  Sequence contains unexpected characters: {sorted(seq_set)}")
            print("Defaulting to amino acid interpretation.")
            return False, SequenceTypeDetector.AMINO_ACIDS
