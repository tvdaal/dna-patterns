"""This module supports analyses of DNA sequences.
"""


import sys
from typing import Optional, List, Tuple, Dict, Any, Callable, Iterable


def parse_txt_file(path: str) -> List[str]:
    """Parse contents of input text file, separated by spaces.

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

    def __init__(self, sequence: str):
        """Initializes the Sequence class.

        Args:
            sequence: Input DNA sequence.
        """
        self.sequence = sequence
        self.length = len(sequence)

    def reverse_complement(self) -> str:
        """Finds the reverse complement of a sequence.

        This method has linear time complexity.

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

    def pattern_count(self, pattern: str) -> int:
        """Counts the number of times a pattern occurs in the sequence.

        This method has linear time complexity.

        Args:
            pattern: The pattern that needs to be matched.

        Returns:
            The number of times the pattern occurs in the sequence.
        """

        count = 0
        n = self.length
        k = len(pattern)
        for i in range(n-k+1):
            if self.sequence[i:i+k] == pattern:
                count += 1

        return count

    def frequency_table(self, k: int) -> Dict[str, int]:
        """Counts the frequency of all k-mers that occur in the sequence.

        This method has linear time complexity.

        Args:
            k: Length of pattern (k-mer).

        Returns:
            A dictionary that contains the frequencies for each k-mer in text.
        """

        frequency_map = {}
        n = self.length
        for i in range(n-k+1):
            pattern = self.sequence[i:i+k]
            count = frequency_map.get(pattern)
            if count is None:
                frequency_map[pattern] = 1
            else:
                count += 1
                frequency_map[pattern] = count

        return frequency_map

    def frequent_words(self, k: int) -> List[str]:
        """Searches for the most frequent patterns of length k in the sequence.

        This method has linear time complexity.

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
