import unittest
from result_process import DirectionHandler, get_coverage, PlotResult


class TestDirectionHandler(unittest.TestCase):
    def setUp(self):
        self.mix_result = {
            ("read1", "ref1", False): [(2, 2)],
            ("read2", "ref2", True): [(3, 1)]  # reverse 3 to forward 8
        }
        self.ref_sequence = {"ref1": "ATCGATCGATCG", "ref2": "GGGCCCAGTACG"}

    def test_sort_direction(self):
        forward, reverse = DirectionHandler.sort_direction(self.mix_result)
        self.assertIn(("read1", "ref1", False), forward)
        self.assertIn(("read2", "ref2", True), reverse)

    def test_merge_direction(self):
        expected_result = {
            ("read1", "ref1", False): [(2, 2)],
            ("read2", "ref2", True): [(8, 1)]  # reverse 3 to forward 8
        }
        merged_result = DirectionHandler.merge_direction(self.mix_result, self.ref_sequence)
        self.assertEqual(dict(merged_result), expected_result)


class TestCoverageCalculation(unittest.TestCase):
    def setUp(self):
        self.unified_result = {
            ("read1", "ref1", False): [(0, 0)],
            ("read2", "ref1", True): [(0, 0)],
            ("read2", "ref2", False): [(5, 0)]
        }
        self.reference_seqs = {"ref1": "ATCGATCGCTCC", "ref2": "GGGCCCGATACC"}
        self.reads_seqs = {"read1": "ATCGAT", "read2": "CGA"}

    def test_coverage_calculation(self):
        coverage = get_coverage(self.unified_result, self.reference_seqs, self.reads_seqs)
        self.assertIn("ref1", coverage)
        self.assertIn("ref2", coverage)
        self.assertEqual(coverage["ref1"], 0.625)  # expected 0.625
        self.assertEqual(coverage["ref2"], 0.125)  # expected 0.125


class TestPlotResult(unittest.TestCase):
    def setUp(self):
        self.match_result = {
            ("read1", "ref1", False): [(100, 2), (1000, 2)],
            ("read2", "ref2", False): [(50, 1), (500, 1), (5000, 1)]
        }

    def test_calculate_match_counts(self):
        match_counts = PlotResult.calculate_match_counts(self.match_result)
        self.assertEqual(match_counts["read1"], 2)
        self.assertEqual(match_counts["read2"], 3)


if __name__ == '__main__':
    unittest.main()
