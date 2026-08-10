import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
import csv


class DirectionHandler:
    """
    To either separate match position to forward strand and reversed strand, or merged them together and
    converts the positions to unified position (forward direction).
    """
    @staticmethod
    def sort_direction(mix_result):
        """
        To separate match position to forward strand and reversed strand with their own 5'-3' positions.
        Arg:
            mix_result: result that contains both forward strand and reversed strand position
        Return:
            forward, reverse: separated results
        """
        forward = defaultdict()
        reverse = defaultdict()
        for (read_id, ref_id, is_reverse), matches in mix_result.items():
            if is_reverse is False:
                forward[(read_id, ref_id, is_reverse)] = matches
            else:
                reverse[(read_id, ref_id, is_reverse)] = matches
        return forward, reverse

    @classmethod
    def merge_direction(cls, mix_result, ref_sequence):
        """
        Convert the positions on reversed strand to it corresponding positions on forward strand,
        and merge it with positions on forward strand.
        Arg:
            mix_result: result that contains both forward strand and reversed strand position
            ref_sequence: reference sequences
        Return:
            forward_unified:
        """
        forward, reverse = cls.sort_direction(mix_result)
        forward_unified = defaultdict(list)

        for (read_id, ref_id, is_reverse), matches in reverse.items():
            ref_length = len(ref_sequence[ref_id])
            for binding_start, hamming_distance in matches:
                unified_start = ref_length-binding_start-1
                forward_unified[(read_id, ref_id, True)].append((unified_start, hamming_distance))
        forward_unified.update(forward)
        return forward_unified


def get_coverage(min_hamming_matcher, reference_seqs, reads_seqs):
    """
    Get reference coverage for each reference sequence
    Arg:
        min_hamming_matcher:
        reference_seqs:
        reads_seqs:
    """
    read_match_sum = defaultdict(int)
    read_to_ref_sum = defaultdict(int)
    coverage = defaultdict(float)

    ref_length = {seq_id: len(seq) for seq_id, seq in reference_seqs.items()}
    read_length = {seq_id: len(seq) for seq_id, seq in reads_seqs.items()}

    for (read_id, ref_id, is_reverse), binding_info in min_hamming_matcher.items():
        read_match_sum[read_id] += len(binding_info)
        read_to_ref_sum[(read_id, ref_id)] += len(binding_info)
    for ref in reference_seqs.keys():
        l_ref = (1 / ref_length[ref])
        for read in reads_seqs.keys():
            if (read, ref) in read_to_ref_sum:
                mij = read_to_ref_sum[(read, ref)]
                mi = read_match_sum[read]
                r = read_length[read]
                coverage[ref] += l_ref * r * mij / mi
    return coverage


class PlotResult:
    def __init__(self, match_result, output_name):
        """
        Arg:
            match_result: assign attribute to self instance
            output_name: prefix for output file
        """
        self.match_counts = self.calculate_match_counts(match_result)
        self.plot_histogram(self.match_counts, output_name)

    @staticmethod
    def calculate_match_counts(match_result):
        """ calculate match positions for each reads """
        match_counts = defaultdict(int)
        for (read_id, ref_id, is_reverse), matches in match_result.items():
            match_counts[read_id] += len(matches)
        return match_counts

    @staticmethod
    def plot_histogram(match_counts, output_name, x_range=None, x_step=None):
        """ Generate the histogram """
        plt.figure(figsize=(10, 6))
        plt.hist(match_counts.values(), bins=40, edgecolor='black', alpha=0.7)
        ax = plt.gca()
        if x_range:
            ax.set_xlim(x_range)
        if x_step:
            ax.xaxis.set_major_locator(ticker.MultipleLocator(x_step))
        plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        plt.xlabel("Nr. best locations")
        plt.ylabel("Count")
        plt.title("Histogram of best location per read")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.savefig(f"{output_name}_best_location_per_read.pdf")
        plt.show(block=False)


class GenerateCsv:
    @staticmethod
    def get_match_csv(best_locations, output_name):
        with open(f"{output_name}_MatchPosition.csv", mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Read_id", "Ref_id", "Position", "Substitution(s)"])
            for (read_id, ref_id, is_reverse), binding_info in best_locations.items():
                for binding_start, hamming_distance in binding_info:
                    writer.writerow([read_id, ref_id, binding_start + 1, hamming_distance])

    @staticmethod
    def get_coverage_csv(ref_coverage, output_name):
        with open(f"{output_name}_ReferenceCoverage.csv", mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Reference", "Coverage"])
            for reference, coverage in ref_coverage.items():
                writer.writerow([reference, coverage])
