from file_handling import FileHandling
from analysis_methods import BuildKmerBinaryIndex, KmerIndexMatcher, HammingDistanceMatcher
from result_process import DirectionHandler
from result_process import PlotResult as GetPlot, GenerateCsv as GetCsv
from result_process import get_coverage
import argparse
import os

"""
*****Commands example and explanation:*****
python3 scripts/main.py --k 15 refs0.fa reads0.fa Refs0_Reads0----Example command
python3 scripts/main.py                                ---interpreter and main function
                --k 15                                 ---length of kmer (default 10)
                       refs1.fa                        ---reference fasta filename
                                reads1.fa              ---read fasta filename
                                          Refs1_Reads1 ---output file prefix
=======================================================================================
!!!!!Important TIPS!!!!!
Filename extensions for the input files are needed
Do NOT use invalid character for the output file prefix! E.g. '/' or  '.'
Length of kmer should be between 5 and 32. Preferably 10~30!
=======================================================================================
*****About INPUT and OUTPUT file*****
Before running the program, place references in 'data/references' and reads in 'data/reads'.
All results will be stored in 'results/output'.

Test code: 
python3 scripts/main.py --k 15 refs0.fa reads0.fa Refs0_Reads0
python3 scripts/main.py --k 15 refs1.fa reads1.fa Refs1_Reads1
python3 scripts/main.py --k 15 refs5.fa reads5.fa Refs5_Reads5
python3 scripts/main.py --k 15 refs10.fa reads10.fa Refs10_Reads10

***To see what will happen when binding locations exceed the ends of reference sequence.
python3 scripts/main.py --k 10 refs5.fa reads5.fa test_test
***To see what will happen if length of k-mer is to large (Error). (For x64 bit processor, k_max = 32)
python3 scripts/main.py refs5.fa reads5.fa test_test
"""


def get_user_input():
    """
    Parses command-line arguments for the program.
    Combine arguments to get data directory, output directory and output absolute path.
    Return:
        reference_path: path to find the reference fasta
        reads_path: path to find the reads fasta
        output_path: prefix of all output files
        k_length: len of the kmer
    """
    parser = argparse.ArgumentParser(description="DNA Short Reads Mapping")
    # mandatory arguments
    parser.add_argument("reference", type=str, help="ref fasta")
    parser.add_argument("reads", type=str, help="reads fasta")
    parser.add_argument("output", type=str, help="file output prefix")
    # optional arguments
    parser.add_argument("--k", type=int, default=10, help="length of k-mer (default 10, max 32)")
    # obtain arguments
    commands = parser.parse_args()
    # set directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(project_dir, "data")
    output_dir = os.path.join(project_dir, "results", "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print('Folder "results/output" has been created to store output files')
    # get complete, absolute path
    reference_path = os.path.join(data_dir, "references", commands.reference)
    reads_path = os.path.join(data_dir, "reads", commands.reads)
    output_path = os.path.join(output_dir, commands.output)
    k_length = commands.k
    return reference_path, reads_path, output_path, k_length


def main():
    # Get user input parameters
    reference_path, reads_path, output_path, k_length = get_user_input()

    #  load fasta file
    reference_seqs = FileHandling(reference_path).load_fasta()
    reads_seqs = FileHandling(reads_path).load_fasta()

    #  Create binary index for sequences
    reference_forward_index = BuildKmerBinaryIndex(reference_seqs, k_length, build_reverse=False).build_binary_index()
    reference_reverse_index = BuildKmerBinaryIndex(reference_seqs, k_length, build_reverse=True).build_binary_index()
    reads_index = BuildKmerBinaryIndex(reads_seqs, k_length, build_reverse=False).build_binary_index()

    #  Matching sequences with the binary index
    kmer_matcher = KmerIndexMatcher(reference_forward_index, reference_reverse_index, reads_index)
    kmer_match_result = kmer_matcher.get_match()
    kmer_match_locations = kmer_matcher.get_locations(kmer_match_result)

    #  Refine matching results with the hamming distance
    hamming_matcher = HammingDistanceMatcher(kmer_match_locations, reference_seqs, reads_seqs)
    hamming_match_locations = hamming_matcher.get_hamming_distances()
    hamming_best_locations = hamming_matcher.get_min_distances(hamming_match_locations)

    #  Merged result and get reference coverage
    unified_result = DirectionHandler.merge_direction(hamming_best_locations, reference_seqs)
    ref_coverage = get_coverage(hamming_best_locations, reference_seqs, reads_seqs)

    # Get the output
    GetPlot(unified_result, output_path)
    GetCsv.get_coverage_csv(ref_coverage, output_path)
    GetCsv.get_match_csv(hamming_best_locations, output_path)
    print("Analysis complete! Check results in folder 'results/output'. ")
    exit(0)


if __name__ == '__main__':
    main()
