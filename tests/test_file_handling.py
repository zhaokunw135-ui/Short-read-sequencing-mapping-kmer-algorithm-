import unittest
import os
from file_handling import FileHandling


class TestFileHandling(unittest.TestCase):

    def setUp(self):
        self.valid_fasta = "valid.fasta"
        with open(self.valid_fasta, 'w') as file:
            file.write(">seq1\n")
            file.write("ATGCAAGAGACGACGTGC\n")
            file.write(">seq2\n")
            file.write("GATTATAGCATCGATACCA\n")

        self.invalid_format_fasta = "invalid_format.fasta"
        with open(self.invalid_format_fasta, 'w') as file:
            file.write("ATGCATGC\n")
            file.write("ATGCATGC\n")

        self.invalid_sequence_fasta = "invalid_seq.fasta"
        with open(self.invalid_sequence_fasta, 'w') as file:
            file.write(">seq1\n")
            file.write("ATGCXADBYZ\n")
            file.write(">seq2\n")
            file.write("GATTACA\n")

        self.no_exist_fasta = "abcdefg.fasta"

    def tearDown(self):
        for file in [self.valid_fasta, self.invalid_format_fasta, self.invalid_sequence_fasta]:
            if os.path.exists(file):
                os.remove(file)

    def test_file_exists(self):
        FileHandling(self.valid_fasta)
        self.assertRaises(FileNotFoundError, FileHandling, self.no_exist_fasta)

    def test_file_format(self):
        FileHandling(self.valid_fasta)
        self.assertRaises(ValueError, FileHandling, self.invalid_format_fasta)

    def test_file_sequence(self):
        FileHandling(self.valid_fasta)
        self.assertRaises(ValueError, FileHandling, self.invalid_sequence_fasta)

    def test_load_fasta(self):
        expected_result = {
            "seq1": "ATGCAAGAGACGACGTGC",
            "seq2": "GATTATAGCATCGATACCA"
        }
        sequences = FileHandling(self.valid_fasta).load_fasta()

        self.assertEqual(sequences, expected_result)


if __name__ == '__main__':
    unittest.main()
