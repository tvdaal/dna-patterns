#!/usr/bin/env python3
"""This module supports analyses of DNA sequences.
"""


import numpy as np
from typing import Optional, List, Tuple, Dict


def parse_txt_file(path: str) -> List[str]:
    """Parses contents of input text file, separated by spaces.

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
    """This class defines a DNA sequence.

    Attributes:
        complements: Dictionary of complementary nucleotides.
    """

    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def __init__(self, sequence: str) -> None:
        """Initializes the Sequence class.

        Args:
            sequence: Input DNA sequence.
        """
        self.sequence = sequence
        self.length = len(sequence)

    def reverse_complement(self) -> str:
        """Finds the reverse complement of a sequence.

        Returns:
            The reverse complement of the sequence.
        """

        rc_list = []
        for nucleotide in self.sequence:
            complement = self.complements[nucleotide]
            rc_list.append(complement)
        rc_list.reverse()
        rc = "".join(rc_list)

        return rc

    def pattern_count(self, pattern: str) -> Tuple[int, List[int]]:
        """Counts the number of times a pattern occurs in the sequence.

        Args:
            pattern: The pattern that needs to be matched.

        Returns:
            The number of times the pattern occurs in the sequence.
        """

        count = 0
        n = self.length
        k = len(pattern)
        positions = []
        for i in range(n-k+1):
            if self.sequence[i:i+k] == pattern:
                count += 1
                positions.append(i)

        return count, positions

    def frequency_table(
        self,
        k: int,
        window: Optional[str] = None,
    ) -> Dict[str, int]:
        """Counts the frequency of all k-mers that occur in the sequence.

        Args:
            k: Length of pattern (k-mer).
            window: A subsequence that can be looked at instead of the full
                sequence.

        Returns:
            A dictionary that contains the frequencies for each k-mer in the
            sequence.
        """

        frequency_map = {}
        if window is None:
            n = self.length
            sequence = self.sequence
        else:
            n = len(window)
            sequence = window

        for i in range(n-k+1):
            pattern = sequence[i:i+k]
            count = frequency_map.get(pattern)
            if count is None:
                frequency_map[pattern] = 1
            else:
                count += 1
                frequency_map[pattern] = count

        return frequency_map

    def frequent_words(self, k: int) -> List[str]:
        """Searches for the most frequent patterns of length k in the sequence.

        Args:
            k: Length of pattern (k-mer).

        Returns:
            A list of most frequent k-mers that appear in the sequence.
        """

        frequent_patterns = []
        frequency_map = self.frequency_table(k)
        max_count = max(frequency_map.items(), key = lambda x: x[1])[1]

        for key, value in frequency_map.items():
            if value == max_count:
                frequent_patterns.append(key)

        return frequent_patterns

    def find_clumps(self, k: int, L: int, t: int) -> List[str]:
        """Searches for clumps of k-mers in the sequence.

        Args:
            k: Length of pattern (k-mer).
            L: Length of window in the sequence that is looked at.
            t: Minimum number of occurrences of the k-mer in the given window.

        Returns:
            A list of distinct k-mers that appear a minimum of t times in any
            of the windows of length L.
        """

        patterns = []
        n = self.length
        for i in range(n-L+1):
            window = self.sequence[i:i+L]
            frequency_map = self.frequency_table(k, window=window)
            for key, _ in frequency_map.items():
                if frequency_map[key] >= t:
                    patterns.append(key)
        patterns = list(set(patterns))

        return patterns

    def skew_graph(self) -> Tuple[List[int], List[int]]:
        """Finds the G-C skew graph for the input sequence, incl. the minima.

        Returns:
            A list of skew scores for each position in the sequence.
            Additionaly, a list of skew minima is returned.
        """

        skew_scores = [0]  # The skew graph starts at 0 by convention.
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

        return skew_scores, skew_minima
