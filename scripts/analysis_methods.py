from collections import defaultdict
from tqdm import tqdm
import numpy as np
import re


class BuildKmerBinaryIndex:
    """
    Generate hash index for sequence
    input format: sequence dict {header:sequence,......}
    output format: index list, [(header, position, binary_kmer),.....]
    """
    def __init__(self, sequence_dict, k_length=10, build_reverse=False):
        """
        Arg:
            sequence_dict: dictionary contains header and sequences
            k_length: the length of the kmer (user input or default value 10)
            build_reverse: only build the reversed reference sequence strand when True
        """
        self.k_length = k_length
        self.sequence_dict = sequence_dict
        self.is_reversed = build_reverse

    def build_binary_index(self, step=1):
        """
        Arg:
            step: step when generating kmers from the sequence, default 1
        Return:
            index: a list, structure (header, position, kmer)
        """
        index = list()

        total_sequences = len(self.sequence_dict)
        progress_bar = tqdm(total=total_sequences, desc=".....Building kmer index")

        for header, sequence in self.sequence_dict.items():
            progress_bar.set_description(f"......Processing: {header}")
            if self.is_reversed:
                sequence = self.to_reverse(sequence)

            for position in range(0, len(sequence) - self.k_length + 1, step):
                kmer = sequence[position:position + self.k_length]
                index.append((header, position, self.to_hash(kmer)))

            progress_bar.update(1)
        progress_bar.set_description("Index building complete")
        progress_bar.close()
        return index

    @staticmethod
    def to_reverse(sequence):
        """
        Calculate the reversed complementary strand
        """
        base_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(base_pair[c] for c in sequence[::-1])

    @staticmethod
    def to_hash(kmer):
        """
        Covert a kmer to its hash value
        """
        base_to_binary = {'A': "00", 'T': "01", 'C': "10", 'G': "11"}
        return int("".join(base_to_binary[base] for base in kmer), 2)


class KmerIndexMatcher:
    """
    Match sequences based on the hash value
    input format:  index list, [(header, position, binary_kmer),.....]
    output format: location dict, {(read_id, ref_id, is_reverse): [binding locations, ......]}
    """
    all_id_maps = {}  # a class variant, stores conversion between numeric id and original id of reference

    def __init__(self, reference_forward_index, reference_reverse_index, reads_index):
        self.reads_index = reads_index
        #  converted reference index to numpy array
        self.reference_forward_index = self.id_to_number(reference_forward_index, "forward")
        self.reference_reverse_index = self.id_to_number(reference_reverse_index, "reverse")

    @classmethod
    def id_to_number(cls, sequence_list, map_key):
        """
        Convert the original id to numeric id in reference sequences and create a contiguous memory
        for reference index, which allow the computer to quickly find the target reference id.
        Arg:
            sequence_list: [(header, position, kmer_hash),...]
            map_key: key to find the id_map
        Return:
            numpy_index: converted sequence_list, id_map
        """
        id_map = {}
        converted = list()
        total_header = len(sequence_list)
        progress_bar = tqdm(total=total_header, desc="......Building mapping for reference", unit="sequence")

        for header, position, kmer_hash in sequence_list:
            numeric_id = int(re.search(r'\d+', header)[0])
            id_map[numeric_id] = header
            converted.append((numeric_id, position, kmer_hash))
            progress_bar.update(1)
        numpy_index = np.array(
            converted,  
            dtype=[('ref_id', 'i4'), ('ref_start', 'i4'), ('ref_kmer_hash', 'u8')]
        )
        cls.all_id_maps[map_key] = id_map
        progress_bar.set_description("ID mapping complete")
        progress_bar.close()
        return numpy_index

    @classmethod
    def number_to_id(cls, numeric_id, map_key):
        """
        Convert the numeric id back to original id in reference sequences
        Arg:
            numeric_id: numeric id of the reference sequence
            map_key: key to find the id_map
        Return:
            the original id of the reference sequence
        """
        id_map = cls.all_id_maps.get(map_key, {})
        return id_map.get(numeric_id)

    def get_match(self):
        """
        Process kmer matching with Numpy vectors
        Return:
            match_result: a dict, (read_id, ref_id, is_reverse) as key, (read_start, ref_start) as value
        """
        match_result = defaultdict(list)

        total_kmers = len(self.reads_index)
        progress_bar = tqdm(total=total_kmers, desc="......K-mer mapping processing ", unit="kmer")

        for read_id, read_start, read_kmer_hash in self.reads_index:
            mask = (self.reference_forward_index['ref_kmer_hash'] == read_kmer_hash)
            matched_refs = self.reference_forward_index[mask]
            mask_rev = (self.reference_reverse_index['ref_kmer_hash'] == read_kmer_hash)
            matched_refs_rev = self.reference_reverse_index[mask_rev]

            for match in matched_refs:
                original_ref_id = self.number_to_id(int(match['ref_id']), "forward")
                match_result[(read_id, original_ref_id, False)].append((read_start, (match['ref_start'])))
            for match in matched_refs_rev:
                original_ref_id = self.number_to_id(int(match['ref_id']), "reverse")
                match_result[(read_id, original_ref_id, True)].append((read_start, (match['ref_start'])))
            progress_bar.update(1)

        progress_bar.set_description("K-mer mapping complete")
        progress_bar.close()

        return match_result

    @staticmethod
    def get_locations(match_result):
        """
        Calculate the theoretical binding start between matched ref and read sequence
        Sort the top_n most frequent binding start locations
        Arg:
            match_result: a dict, (read_id, ref_id, is_reverse) as key, (read_start, ref_start) as value
            locations:  a dict, take (read_id, ref_id,  is_reverse) as key and a list of binding locations as value
        """
        match_locations = defaultdict(set)

        total_reads = len(match_result)
        progress_bar = tqdm(total=total_reads, desc="......Computing locations", unit="read")

        for seq_info, positions in match_result.items():
            for read_start, ref_start in positions:
                binding_start = ref_start - read_start
                match_locations[seq_info].add(int(binding_start))

            progress_bar.update(1)

        progress_bar.set_description("Locations obtained")
        progress_bar.close()

        return match_locations


class HammingDistanceMatcher:
    """
    Get hamming distance between each matched ref and read sequence.
    input format: location dict, {(read_id, ref_id, is_reverse): [binding locations,......]}
    output format: location dict, {(read_id, ref_id, is_reverse): (binding location，hamming distance),......}
    """
    def __init__(self, match_locations, reference_seqs, reads_seqs):
        self.match_locations = match_locations
        self.reference_seqs = reference_seqs
        self.reads_seqs = reads_seqs

    @staticmethod
    def get_reference_fragment(ref_id, ref_seq, read_id, read_seq, start_pos):
        """
        Segment the reference sequence to a fraction that has the same length with the corresponding read sequence,
        based on the start position and the length of read. If the start position locates outside the reference
        sequence, add X letter to the end(or start) of fraction to get the same length with read.
        Arg:
            ref_id: reference sequence header
            ref_seq: reference sequence
            read_id: read sequence header
            read_seq: read sequence
            start_pos: start position
            progress_bar: for displaying the progress
        """
        read_length = len(read_seq)
        ref_length = len(ref_seq)

        end_pos = start_pos + read_length
        if start_pos < 0:
            ref_fragment = 'X' * abs(start_pos) + ref_seq[:end_pos]
            print(f"Binding start: '{start_pos+1}' from '{read_id}' to '{ref_id}' exceed the 5' start of the reference")
        elif end_pos > ref_length:
            ref_fragment = ref_seq[start_pos:] + 'X' * (end_pos - ref_length)
            print(f"Binding start: '{end_pos+1}' from '{read_id}' to  '{ref_id}' exceed the 3' end of the reference")
        else:
            ref_fragment = ref_seq[start_pos:end_pos]

        return ref_fragment

    def get_hamming_distances(self):
        """
        Calculate Hamming distances
        """
        hamming_distances = defaultdict(list)
        total_pairs = len(self.match_locations)
        progress_bar = tqdm(total=total_pairs, desc="......Calculating Hamming distance", unit="comparison")

        for (read_id, ref_id, is_reverse), binding_starts in self.match_locations.items():
            read_seq = self.reads_seqs[read_id]
            ref_seq = self.reference_seqs[ref_id]
            for binding_start in binding_starts:
                ref_frag = self.get_reference_fragment(ref_id, ref_seq, read_id, read_seq, binding_start)
                hamming_distance = sum(base1 != base2 for base1, base2 in zip(ref_frag, read_seq))
                hamming_distances[(read_id, ref_id, is_reverse)].append((binding_start, hamming_distance))
                progress_bar.update(1)

        progress_bar.set_description("Hamming distance calculation complete")
        progress_bar.close()
        return hamming_distances

    @staticmethod
    def get_min_distances(hamming_distances, min_distance=5):
        """
        Filter in the distances no longer than min_distance
        Arg:
            hamming_distances: unsorted result
            min_distance: the minimum allowed distance for sorting
        Return:
            min_hamming_distances: dict, {(read_id, ref_id, is_reverse): hamming distance,......}
        """
        min_hamming_distances = defaultdict(list)
        for id_info, binding_info in hamming_distances.items():
            for binding_start, hamming_distance in binding_info:
                if hamming_distance <= min_distance:
                    min_hamming_distances[id_info].append((binding_start, hamming_distance))

        return min_hamming_distances
