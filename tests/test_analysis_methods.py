import unittest
from analysis_methods import BuildKmerBinaryIndex, KmerIndexMatcher, HammingDistanceMatcher


class TestAllMatching(unittest.TestCase):
    def setUp(self):
        self.reference_seqs = {
            "ref1": "ATCGATCGATCGA",
        }
        self.reads_seqs = {
            "read1": "ATCGATCGATCA",
            "read2": "CGATCGATCGA",
        }

        self.reference_index = BuildKmerBinaryIndex(self.reference_seqs, k_length=11).build_binary_index()
        self.reads_index = BuildKmerBinaryIndex(self.reads_seqs, k_length=11).build_binary_index()

        self.matcher = KmerIndexMatcher(self.reference_index, [], self.reads_index)
        self.match_result = self.matcher.get_match()

        self.match_locations = self.matcher.get_locations(self.match_result)

    def test_to_reverse(self):
        original_seq = "ATGCATGCATGCATGC"
        expected_result = "GCATGCATGCATGCAT"
        self.assertEqual(BuildKmerBinaryIndex.to_reverse(original_seq), expected_result)

    def test_to_hash(self):
        original_seq = "AAAAAAAAAAAG"   # A = 00, G = 11
        self.assertEqual(BuildKmerBinaryIndex.to_hash(original_seq), 3)

    def test_build_kmer_index(self):
        expected_result = [
            ('ref1', 0, 444102),  # k-mer = ATCGATCGATC
            ('ref1', 1, 1776411),  # k-mer = TCGATCGATCG
            ('ref1', 2, 2911340)  # k-mer = CGATCGATCGA
        ]

        expected_result_read = [
            ('read1', 0, 444102),  # k-mer = ATCGATCGATC
            ('read1', 1, 1776408),  # k-mer = TCGATCGATCA
            ('read2', 0, 2911340)   # k-mer = CGATCGATCGA
        ]
        self.assertEqual(self.reference_index, expected_result)
        self.assertEqual(self.reads_index, expected_result_read)

    def test_get_match(self):
        expected_result = {
            ('read1', "ref1", False): [(0, 0)],  # k-mer = ATCGATCGATC
            ('read2', "ref1", False): [(0, 2)]  # k-mer = CGATCGATCGA
        }
        self.assertEqual(self.match_result, expected_result)

    def test_hamming_distance_calculation(self):
        expected_result = {
            ('read1', "ref1", False): [(0, 1)],
            ('read2', "ref1", False): [(2, 0)]
        }
        hamming_matcher = HammingDistanceMatcher(self.match_locations, self.reference_seqs, self.reads_seqs)
        hamming_distances = dict(hamming_matcher.get_hamming_distances())
        for position in hamming_distances.values():
            for read_location, ref_location in position:
                int(read_location)
                int(ref_location)

        self.assertEqual(hamming_distances, expected_result)

    def test_min_hamming_distances(self):
        hamming_matcher = HammingDistanceMatcher(self.match_locations, self.reference_seqs, self.reads_seqs)
        hamming_distances = hamming_matcher.get_hamming_distances()
        min_distances = hamming_matcher.get_min_distances(hamming_distances, min_distance=1)

        for key, values in min_distances.items():
            for seqs, distance in values:
                self.assertLessEqual(distance, 1)


if __name__ == '__main__':
    unittest.main()
