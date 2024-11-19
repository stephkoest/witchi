import numpy as np
from Bio import AlignIO

class AlignmentReader:
    def __init__(self, file, format='fasta'):
        self.file = file
        self.format = format

    def run(self):
        """Main method to read and detect alignment."""
        alignment = self.read_alignment()
        alignment_array = self.make_alignment_array(alignment)
        return alignment, alignment_array

    def read_alignment(self):
        """Read the alignment file."""
        return AlignIO.read(self.file, self.format)

    def make_alignment_array(self, alignment):
        """Convert the alignment to a numpy array."""
        alignment_array = np.array([list(record.seq) for record in alignment])
        alignment_array = np.char.upper(alignment_array)
        return alignment_array