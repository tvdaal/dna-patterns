#!/usr/bin/env python3
"""This module supports analyses of DNA sequences.

Its main object is the Sequence class. This class contains many methods to
study a given DNA sequence.

This module cannot be run directly, but should be imported to another module
or script.

  Typical usage example for employing the Sequence class in a different script:

  import sequence as seq
  import sys
  input_path = sys.argv[1]
  contents = seq.parse_fasta_file(input_path)
  dna_sequence = seq.Sequence(contents[0])
  rev_complement = dna_sequence.reverse_complement()
"""


from __future__ import annotations  # Needed for Python 3.10 type annotations for classes
import itertools
import numpy as np
import random
import sys
from typing import Optional, List, Tuple, Dict, Set, Union, Any


nucleotides = {"A", "T", "G", "C"}
complements = {"A": "T", "T": "A", "G": "C", "C": "G"}


def parse_txt_file(path: str) -> List[str]:
    """Parses contents of an input text file, separated by spaces.

    Args:
        path: Path to input text file.

    Returns:
        The contents of the text file as a list of strings.
    """

    file = open(path, "r")
    contents = file.read()
    text = contents.split()
    file.close()

    return text


def parse_fasta_file(path: str) -> str:
    """Parses contents of an input text file in FASTA format.

    Args:
        path: Path to input text file.

    Returns:
        The sequence.
    """

    file = open(path, "r")
    contents = file.readlines()
    text = "".join(contents[1:]).replace("\n", "")
    file.close()

    return text


class Sequence:
    """This class defines a sequence of nucleotides, i.e. a piece of DNA.

    The meaning of the sequence can range from a simple pattern to an entire
    genome. The sequence is assumed to be of string format.
    """

    def __init__(self, sequence: str) -> None:
        """Initializes the Sequence class.

        Args:
            sequence: Input sequence.
        """
        self.sequence = sequence
        self.length = len(sequence)

    def reverse_complement(self) -> Sequence:
        """Finds the reverse complement of the sequence.

        Returns:
            The reverse complement of the sequence.
        """

        rev_comp_list = []
        for nucleotide in self.sequence:
            complement = complements[nucleotide]
            rev_comp_list.append(complement)
        rev_comp_list.reverse()
        rev_comp = "".join(rev_comp_list)
        rev_comp = Sequence(rev_comp)

        return rev_comp

    def pattern_count(
        self,
        pattern: Sequence,
        hamming_max: int = 0,
    ) -> Dict[str, Union[int, List[int]]]:
        """Counts the number of times a given pattern occurs in the sequence.

        Additionally, the positions of the matches are provided. Overlapping
        patterns are allowed, as well as imperfect matches (controlled by the
        hamminx_max parameter).

        Args:
            pattern: The pattern that needs to be matched.
            hamming_max: The maximum number of allowed mismatches with respect
                to the given pattern.

        Returns:
            The number of times the pattern occurs in the sequence, as well as
            the positions of where they occur.
        """

        seq = self.sequence
        pat = pattern.sequence

        count = 0
        positions = []
        last_pos = len(seq) - len(pat) + 1
        for i in range(last_pos):
            sub_seq = seq[i:i+len(pat)]
            sub_seq = Sequence(sub_seq)
            hamming_dist = sub_seq.hamming_distance(pattern)
            if hamming_dist <= hamming_max:
                count += 1
                positions.append(i)
        results = {"Count": count, "Positions": positions}

        return results

    def hamming_distance(self, sequence: Sequence) -> int:
        """Counts the number of mismatches between two equal-length sequences.

        Args:
            sequence: The sequence against which a comparison is made.

        Returns:
            The number of mismatches between the two sequences, also referred
            to as the Hamming distance.

        Raises:
            ValueError: If the sequences are not of equal length.
        """

        seq_1 = self.sequence
        seq_2 = sequence.sequence
        if len(seq_1) != len(seq_2):
            raise ValueError("The sequences are not of equal length.")

        hamming_distance = 0
        for i in range(len(seq_1)):
            if seq_1[i] != seq_2[i]:
                hamming_distance += 1

        return hamming_distance

    def neighbors(self, hamming_max: int) -> Set[str]:
        """Defines the neighborhood for the sequence.

        A neighborhood consists of all sequences that differ at most by
        hamming_max from the given sequence.

        Args:
            hamming_max: The maximum number of allowed mismatches with respect
                to the sequence.

        Returns:
            A set of patterns that are neighbors to the input sequence.
        """

        seq = self.sequence
        if hamming_max == 0:
            return {seq}
        if self.length == 1:
            return nucleotides

        suffix = seq[1:]
        suffix = Sequence(suffix)
        suffix_neighborhood = suffix.neighbors(hamming_max)

        neighborhood = set()
        for pattern in suffix_neighborhood:
            pattern = Sequence(pattern)
            if suffix.hamming_distance(pattern) < hamming_max:
                for nucleotide in nucleotides:
                    neighborhood.add(nucleotide + pattern.sequence)
            else:
                neighborhood.add(seq[0] + pattern.sequence)

        return neighborhood

    def list_patterns(
        self,
        pattern_length: int,
        unique: bool = True,
    ) -> List[str]:
        """Finds all patterns of given length in a given sequence.

        Args:
            pattern_length: Fixed length of patterns.
            unique: If True only unique patterns are collected.

        Returns:
            List of patterns with length pattern_length.
        """

        patterns = []
        last_pos = self.length - pattern_length + 1
        for i in range(last_pos):
            pattern = self.sequence[i:i+pattern_length]
            if unique:
                if pattern not in patterns:
                    patterns.append(pattern)
            else:
                patterns.append(pattern)

        return patterns

    def frequency_table(
        self,
        pattern_length: int,
        hamming_max: int = 0,
        reverse: bool = True,
    ) -> Dict[str, int]:
        """Counts the frequency of all patterns that occur in the sequence.

        Approximate occurrences of a given pattern can be accounted for too
        with the hamming_max parameter.

        Args:
            pattern_length: Fixed length of patterns that will be looked at.
            hamming_max: The maximum number of allowed mismatches between
                similar patterns.
            reverse: If True then reverse complementary patterns are taken into
                account as well.

        Returns:
            A dictionary that contains the frequencies for each pattern that
            occurs in the sequence.
        """

        frequency_map = {}
        patterns = self.list_patterns(pattern_length)
        for pattern in patterns:
            pattern = Sequence(pattern)
            all_patterns = []
            all_patterns.append(pattern)
            if reverse:
                pattern_rc = pattern.reverse_complement()
                all_patterns.append(pattern_rc)
            for pat in all_patterns:
                neighborhood = pat.neighbors(hamming_max)
                for neighbor in neighborhood:
                    count = frequency_map.get(neighbor)
                    if count is None:
                        frequency_map[neighbor] = 1
                    else:
                        count += 1
                        frequency_map[neighbor] = count

        return frequency_map

    def frequent_words(
        self,
        pattern_length: int,
        hamming_max: int = 0,
        reverse: bool = True,
    ) -> List[str]:
        """Searches for the most frequent patterns in the sequence.

        Approximate occurrences of a given pattern can be accounted for too
        with the hamming_max parameter.

        Args:
            pattern_length: Fixed length of patterns.
            hamming_max: The maximum number of allowed mismatches between
                similar patterns.
            reverse: If True then reverse complementary patterns are taken into
                account as well.

        Returns:
            A list of most frequent patterns of fixed length that appear in the
            sequence.
        """

        frequent_patterns = []
        frequency_map = self.frequency_table(
            pattern_length,
            hamming_max=hamming_max,
            reverse=reverse,
        )
        max_count = max(frequency_map.items(), key = lambda x: x[1])[1]

        for key, value in frequency_map.items():
            if value == max_count:
                frequent_patterns.append(key)

        return frequent_patterns

    def find_clumps(
        self,
        pattern_length: int,
        window_length: int,
        min_occur: int,
        hamming_max: int = 0,
        reverse: bool = True,
    ) -> Dict[str, Union[int, List[str]]]:
        """Searches for clumps of patterns in the sequence.

        Approximate occurrences of a given pattern can be accounted for too
        with the hamming_max parameter.

        Args:
            pattern_length: Fixed length of patterns.
            window_length: Length of window in the sequence that is looked at.
            min_occur: Minimum number of occurrences of the pattern in the
                given window.
            hamming_max: The maximum number of allowed mismatches between
                similar patterns.
            reverse: If True then reverse complementary patterns are taken into
                account as well.

        Returns:
            A list of distinct patterns that appear a minimum of min_occur
            times in any of the windows of given length. Also the number of
            different patterns is given.
        """

        patterns = []
        last_pos = self.length - window_length + 1
        for i in range(last_pos):
            window = self.sequence[i:i+window_length]
            window = Sequence(window)
            frequency_map = window.frequency_table(
                pattern_length,
                hamming_max=hamming_max,
                reverse=reverse,
            )
            for key, _ in frequency_map.items():
                if frequency_map[key] >= min_occur:
                    patterns.append(key)

        patterns = list(set(patterns))
        number_patterns = len(patterns)
        results = {"Number": number_patterns, "Patterns": patterns}

        return results

    def skew_graph(self) -> Dict[str, List[int]]:
        """Finds the G-C skew graph for the input sequence, incl. its minima.

        Returns:
            A list of skew scores for each position in the sequence.
            Additionaly, a list of skew minima is returned.
        """

        skew_scores = [0]  # The skew graph starts with a zero by convention.
        skew_minima = []
        for i in range(self.length):
            if self.sequence[i] == "G":
                score = skew_scores[i] + 1
            elif self.sequence[i] == "C":
                score = skew_scores[i] - 1
            else:
                score = skew_scores[i]
            skew_scores.append(score)

        skew_array = np.array(skew_scores)
        skew_minima = list(np.where(skew_array == skew_array.min())[0])
        results = {"Skew scores": skew_scores, "Skew minima": skew_minima}

        return results

    def most_probable_string(
        self,
        pattern_length: int,
        profile_matrix: Dict[str, List[float]],
    ) -> str:
        """Finds the profile-most probable string in a given sequence.

        Args:
            pattern_length: Fixed length of patterns.
            profile_matrix: The profile matrix for an unspecified collection of
                sequences.

        Returns:
            The most probable string of length pattern_length given the profile
            matrix.
        """

        highest_prob = 0.0
        most_prob_string = self.sequence[:pattern_length]
        patterns = self.list_patterns(pattern_length)
        for pattern in patterns:
            total_prob = 1.0
            for i, nucleotide in enumerate(pattern):
                prob = profile_matrix[nucleotide][i]
                total_prob *= prob
            if total_prob > highest_prob:
                highest_prob = total_prob
                most_prob_string = pattern

        return most_prob_string

    def biased_string(
        self,
        pattern_length: int,
        profile_matrix: Dict[str, List[float]],
    ) -> str:
        """Finds a profile-randomly generated string in a given sequence.

        Args:
            pattern_length: Fixed length of patterns.
            profile_matrix: The profile matrix for an unspecified collection of
                sequences.

        Returns:
            A string of length pattern_length based on a probability
            distribution that derives from the given profile matrix.
        """

        patterns = self.list_patterns(pattern_length, unique=False)
        probs = []
        for pattern in patterns:
            total_prob = 1.0
            for i, nucleotide in enumerate(pattern):
                prob = profile_matrix[nucleotide][i]
                total_prob *= prob
            probs.append(total_prob)

        # Normalize the probability distribution to 1:
        probs_arr = np.array(probs)
        norm_arr = probs_arr / np.sum(probs_arr)
        norm_probs = norm_arr.tolist()

        # Pick a random pattern based on the defined probability distribution:
        biased_string = random.choices(patterns, weights=norm_probs)[0]

        return biased_string


def find_motifs(
    sequences: List[Sequence],
    pattern_length: int,
    hamming_max: int,
) -> Set[str]:
    """Finds regulatory motifs across the sequences in a brute force way.

    The motif has to occur in every single sequence. Minor mismatches are
    allowed and controlled by the hamming_max parameter.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern_length: Fixed length of patterns.
        hamming_max: The maximum number of allowed mismatches between
            similar patterns.

    Returns:
        A set of regulatory motifs.
    """

    motifs = set()
    first_seq = sequences[0]
    last_pos = first_seq.length - pattern_length + 1
    for i in range(last_pos):
        pat = Sequence(first_seq.sequence[i:i+pattern_length])
        neighborhood = pat.neighbors(hamming_max)
        for neighbor in neighborhood:
            pattern = Sequence(neighbor)
            counts = []
            for seq in sequences[1:]:
                count = seq.pattern_count(pattern, hamming_max)["Count"]
                if count == 0:
                    break
                counts.append(count)
            if len(counts) == len(sequences[1:]):  # The pattern then appears in all sequences.
                motifs.add(neighbor)

    return motifs


def distance_pattern_strings(
    sequences: List[Sequence],
    pattern: Sequence,
) -> Dict[str, Union[List[Optional[Sequence]], int]]:
    """Calculates the sum of distances between the pattern and each string.

    The Hamming distance is minimized for each sequence, respecting a fixed
    pattern length.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern: The reference pattern.

    Returns:
        The motif matrix and the total (minimum) Hamming distance between the
        pattern and the collection of sequences (the mofif score).
    """

    best_motifs = []
    total_distance = 0
    for seq in sequences:
        motif = None
        distance = sys.maxsize
        patterns = seq.list_patterns(pattern.length)
        for key in patterns:
            sub_seq = Sequence(key)
            hamming_distance = pattern.hamming_distance(sub_seq)
            if distance > hamming_distance:
                distance = hamming_distance
                motif = sub_seq
        best_motifs.append(motif)
        total_distance += distance
    results = {"Motif matrix": best_motifs, "Score": total_distance}

    return results


def median_motif_search(
    sequences: List[Sequence],
    pattern_length: int,
) -> Dict[str, Union[List[Sequence], Sequence, int]]:
    """Finds the median string for the given collection of sequences.

    The Hamming distance is minimized for each sequence, respecting a fixed
    length of pattern_length.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern_length: Fixed length of patterns.

    Returns:
        The median string (consensus motif) for the collection of sequences,
        as well as the motif matrix and the motif score.
    """

    min_distance = sys.maxsize
    generated_patterns = itertools.product(
        list(nucleotides),
        repeat=pattern_length,
    )
    patterns = ["".join(list(tup)) for tup in list(generated_patterns)]
    for pattern in patterns:
        pattern = Sequence(pattern)
        results = distance_pattern_strings(sequences, pattern)
        total_distance = results["Score"]
        if min_distance > total_distance:
            min_distance = total_distance
            best_motifs = results["Motif matrix"]
            median_string = pattern

    results = {
        "Motif matrix": best_motifs,
        "Score": min_distance,
        "Consensus": median_string,
    }

    return results


def find_consensus_motif(
    motifs: List[Optional[Sequence]],
    motif_length: int,
    laplace: bool = True,
) -> Dict[str, Union[Dict[str, List[float]], int, str]]:
    """Constructs the profile matrix and finds the consensus motif.

    The motifs score is obtained by summing the elements in the motifs
    matrix that do not coincide with the consensus elements.

    Args:
        motifs: Collection of motifs for which a profile matrix will be
            constructed.
        motif_length: Fixed length of motifs.
        laplace: If True then pseudocounts are added to the profile matrix
            to reduce its sparsity.

    Returns:
        A profile matrix for the given collection of motifs, as well as
        the overall motifs score and the consensus motif.
    """

    profile_matrix = {"A": [], "C": [], "G": [], "T": []}
    motifs_score = 0
    consensus_motif = []
    for i in range(motif_length):
        slicing = []
        num_motifs = 0
        for motif in motifs:
            if motif is not None:
                slicing.append(motif.sequence[i])
                num_motifs += 1
        slice_str = "".join(slicing)
        count_max = 0
        score = sys.maxsize
        consensus_nucleotide = None

        for nucleotide in nucleotides:
            count = slice_str.count(nucleotide)
            if laplace:
                frac = (count + 1) / (num_motifs + len(nucleotides))
            else:
                frac = count / num_motifs
            profile_matrix[nucleotide].append(frac)
            if count > count_max:
                count_max = count
                score = num_motifs - count_max
                consensus_nucleotide = nucleotide
        motifs_score += score
        consensus_motif.append(consensus_nucleotide)

    consensus_motif = "".join(consensus_motif)
    results = {
        "Matrix": profile_matrix,
        "Score": motifs_score,
        "Consensus": consensus_motif,
    }

    return results


def greedy_motif_search(
    sequences: List[Sequence],
    pattern_length: int,
    laplace: bool = True,
) -> Dict[str, Union[List[Sequence], Sequence, int]]:
    """Finds a motif matrix in a greedy way.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern_length: Fixed length of patterns.
        laplace: If True then pseudocounts are added to the profile matrix
            to reduce its sparsity.

    Returns:
        A motif matrix for the given collection of sequences, as well as
        the consensus string and the motif score.
    """

    best_motifs = [Sequence(seq.sequence[:pattern_length]) for seq in sequences]
    best_motifs_score = find_consensus_motif(
        best_motifs,
        pattern_length,
        laplace=laplace,
    )["Score"]
    consensus_motif = best_motifs[0]

    num_seqs = len(sequences)
    patterns = sequences[0].list_patterns(pattern_length)  # Only collect patterns from the first sequence.
    for motif in patterns:
        motifs = [None] * num_seqs
        motifs[0] = Sequence(motif)
        for i in range(1, num_seqs):
            profile_matrix = find_consensus_motif(
                motifs,
                pattern_length,
                laplace=laplace,
            )["Matrix"]
            motif_str = sequences[i].most_probable_string(
                pattern_length,
                profile_matrix,
            )
            motifs[i] = Sequence(motif_str)

        motif_results = find_consensus_motif(
            motifs,
            pattern_length,
            laplace=laplace,
        )
        motifs_score = motif_results["Score"]
        if motifs_score < best_motifs_score:
            best_motifs = motifs
            best_motifs_score = motifs_score
            consensus_motif = Sequence(motif_results["Consensus"])

    results = {
        "Motif matrix": best_motifs,
        "Score": best_motifs_score,
        "Consensus": consensus_motif,
    }

    return results


def random_motif_search(
    sequences: List[Sequence],
    pattern_length: int,
    laplace: bool = True,
) -> Dict[str, Union[List[Sequence], Sequence, int]]:
    """Finds a motif matrix in a randomized search.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern_length: Fixed length of patterns.
        laplace: If True then pseudocounts are added to the profile matrix
            to reduce its sparsity.

    Returns:
        A motif matrix for the given collection of sequences, as well as
        the consensus string and the motif score.
    """

    num_seqs = len(sequences)
    last_pos = sequences[0].length - pattern_length + 1
    motifs = [None] * num_seqs
    for i, seq in enumerate(sequences):
        start_pos = random.randrange(0, last_pos)
        end_pos = start_pos + pattern_length
        motifs[i] = Sequence(seq.sequence[start_pos:end_pos])
    best_motifs = motifs
    best_motifs_results = find_consensus_motif(
        best_motifs,
        pattern_length,
        laplace=laplace,
    )
    best_motifs_score = best_motifs_results["Score"]
    consensus_motif = best_motifs_results["Consensus"]

    while True:
        profile_matrix = find_consensus_motif(
            motifs,
            pattern_length,
            laplace=laplace,
        )["Matrix"]
        for i, seq in enumerate(sequences):
            motif_str = seq.most_probable_string(pattern_length, profile_matrix)
            motifs[i] = Sequence(motif_str)

        motif_results = find_consensus_motif(
            motifs,
            pattern_length,
            laplace=laplace,
        )
        motifs_score = motif_results["Score"]
        if motifs_score < best_motifs_score:
            best_motifs = motifs
            best_motifs_score = motifs_score
            consensus_motif = Sequence(motif_results["Consensus"])
        else:
            results = {
                "Motif matrix": best_motifs,
                "Score": best_motifs_score,
                "Consensus": consensus_motif,
            }
            return results


def gibbs_motif_search(
    sequences: List[Sequence],
    pattern_length: int,
    inner_loop_iterations: int,
    laplace: bool = True,
) -> Dict[str, Union[List[Sequence], Sequence, int]]:
    """Finds a motif matrix using Gibbs sampling.

    Args:
        sequences: Collection of sequences that will be searched for
            patterns.
        pattern_length: Fixed length of patterns.
        inner_loop_iterations: Number of times Gibbs sampling takes place.
        laplace: If True then pseudocounts are added to the profile matrix
            to reduce its sparsity.

    Returns:
        A motif matrix for the given collection of sequences, as well as
        the consensus string and the motif score.
    """

    num_seqs = len(sequences)
    last_pos = sequences[0].length - pattern_length + 1
    motifs = [None] * num_seqs
    for i, seq in enumerate(sequences):
        start_pos = random.randrange(0, last_pos)
        end_pos = start_pos + pattern_length
        motifs[i] = Sequence(seq.sequence[start_pos:end_pos])
    best_motifs = motifs
    best_motifs_results = find_consensus_motif(
        best_motifs,
        pattern_length,
        laplace=laplace,
    )
    best_motifs_score = best_motifs_results["Score"]
    consensus_motif = best_motifs_results["Consensus"]

    for j in range(inner_loop_iterations):
        seq_idx = random.randrange(0, num_seqs)
        motifs[seq_idx] = None
        profile_matrix = find_consensus_motif(
            motifs,
            pattern_length,
            laplace=laplace,
        )["Matrix"]
        motif_str = sequences[seq_idx].biased_string(
            pattern_length,
            profile_matrix,
        )
        motifs[seq_idx] = Sequence(motif_str)

        motif_results = find_consensus_motif(
            motifs,
            pattern_length,
            laplace=laplace,
        )
        motifs_score = motif_results["Score"]
        if motifs_score < best_motifs_score:
            best_motifs = motifs
            best_motifs_score = motifs_score
            consensus_motif = Sequence(motif_results["Consensus"])

    results = {
        "Motif matrix": best_motifs,
        "Score": best_motifs_score,
        "Consensus": consensus_motif,
    }

    return results


def repeated_motif_search(
    sequences: List[Sequence],
    pattern_length: int,
    num_iterations: int,
    gibbs: bool = True,
    gibbs_inner_loop: int = 1000,
    laplace: bool = True,
) -> Dict[str, Union[List[Sequence], Sequence, int]]:
    """Finds the best motifs after repeated runs of a motif search algorithm.

    The motif search is repeated num_iterations times upon which the best
    motifs are selected.

    Args:
        sequences: Collection of sequences that will be searched for patterns.
        pattern_length: Fixed length of patterns.
        num_iterations: Number of runs of the chosen algorithm.
        gibbs: If True then the Gibbs sampling algorithm is used. Otherwise the
            randomized search algorithm is used.
        gibbs_inner_loop: Number of times Gibbs sampling takes place.
        laplace: If True then pseudocounts are added to the profile matrix
            to reduce its sparsity.

    Returns:
        A motif matrix for the given collection of sequences, as well as
        the consensus string and the motif score.
    """

    results = {}
    if gibbs:
        for i in range(1, num_iterations+1):
            iter_result = gibbs_motif_search(
                sequences,
                pattern_length,
                gibbs_inner_loop,
                laplace=laplace,
            )
            results.update({i: iter_result})
    else:
        for i in range(1, num_iterations+1):
            iter_result = random_motif_search(
                sequences,
                pattern_length,
                laplace=laplace,
            )
            results.update({i: iter_result})

    best_iter = min(results, key=lambda key: results[key]["Score"])
    best_result = results[best_iter]

    return best_result
