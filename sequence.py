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
  contents = seq.parse_txt_file(input_path)
  sequence = seq.Sequence(contents[0])
  rev_complement = sequence.reverse_complement()
"""


from __future__ import annotations  # Needed for Python 3.10 type annotations for classes
import numpy as np
from typing import Optional, List, Tuple, Dict, Union


def parse_txt_file(path: str) -> List[str]:
    """Parses contents of an input text file, separated by spaces.

    Args:
        path: Path to input text file.

    Returns:
        The contents of the text file as a list of strings.
    """

    f = open(path, "r")
    contents = str(f.read())
    f.close()

    return contents.split()


class Sequence:
    """This class defines a sequence of nucleotides, i.e. a piece of DNA.

    The meaning of the sequence can range from a simple pattern to an entire
    genome. The sequence is assumed to be of string format.

    Attributes:
        complements: Dictionary of complementary nucleotides.
    """

    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def __init__(self, sequence: str) -> None:
        """Initializes the Sequence class.

        Args:
            sequence: Input sequence.
        """
        self.sequence = sequence
        self.length = len(sequence)

    def reverse_complement(self) -> str:
        """Finds the reverse complement of the sequence.

        Returns:
            The reverse complement of the sequence.
        """

        rev_comp_list = []
        for nucleotide in self.sequence:
            complement = self.complements[nucleotide]
            rev_comp_list.append(complement)
        rev_comp_list.reverse()
        rev_comp = "".join(rev_comp_list)

        return rev_comp

    def pattern_count(
        self,
        pattern: Sequence,
        hamming_thr: int = 0,
    ) -> Dict[str, Union[int, List[int]]]:
        """Counts the number of times a given pattern occurs in the sequence.

        Additionally, the positions of the matches are provided. Overlapping
        patterns are allowed, as well as imperfect matches.

        Args:
            pattern: The pattern that needs to be matched.
            hamming_thr: If specified, this threshold provides the maximum
                number of allowed mismatches with respect to the given pattern.

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
            if hamming_dist <= hamming_thr:
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
        """

        seq_1 = self.sequence
        seq_2 = sequence.sequence
        assert len(seq_1) == len(seq_2), "The sequences are not of equal length."

        hamming_distance = 0
        for i in range(len(seq_1)):
            if seq_1[i] != seq_2[i]:
                hamming_distance += 1

        return hamming_distance

    def frequency_table(
        self,
        pattern_length: int,
        window: Optional[Sequence] = None,
    ) -> Dict[str, int]:
        """Counts the frequency of all patterns that occur in the sequence.

        Args:
            pattern_length: Fixed length of patterns that will be looked at.
            window: A subsequence that can be looked at instead of the full
                sequence at hand.

        Returns:
            A dictionary that contains the frequencies for each pattern that
            occurs in the sequence.
        """

        if window is None:
            seq = self.sequence
        else:
            seq = window.sequence
        string_length = len(seq)

        frequency_map = {}
        last_pos = string_length - pattern_length + 1
        for i in range(last_pos):
            pat = seq[i:i+pattern_length]
            count = frequency_map.get(pat)
            if count is None:
                frequency_map[pat] = 1
            else:
                count += 1
                frequency_map[pat] = count

        return frequency_map

    def frequent_words(self, pattern_length: int) -> List[str]:
        """Searches for the most frequent patterns in the sequence.

        Args:
            pattern_length: Fixed length of patterns.

        Returns:
            A list of most frequent patterns of fixed length that appear in the
            sequence.
        """

        frequent_patterns = []
        frequency_map = self.frequency_table(pattern_length)
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
    ) -> Dict[str, Union[int, List[str]]]:
        """Searches for clumps of identical patterns in the sequence.

        Args:
            pattern_length: Fixed length of patterns.
            window_length: Length of window in the sequence that is looked at.
            min_occur: Minimum number of occurrences of the pattern in the
                given window.

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
            frequency_map = self.frequency_table(pattern_length, window=window)
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

        skew_scores = [0]  # The skew graph with a zero by convention.
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
