"""This module supports pattern counting in DNA sequences.

  Example of how to run:

  python pattern_count_1.2.5.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Dataset_1.2.5.txt
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

    return contents.split()


def pattern_count(text: str, pattern: str) -> int:
    """Counts the number of times a pattern occurs in a given string.

    Args:
        text: Input DNA sequence.
        pattern: DNA pattern to match.

    Returns:
        The number of times the pattern occurs in the DNA sequence.
    """

    count = 0
    len_text = len(text)
    len_pattern = len(pattern)
    max_idx = len_text - len_pattern
    for i in range(max_idx+1):
        if text[i:i+len_pattern] == pattern:
            count += 1

    return count


if __name__ == "__main__":
    input_path = sys.argv[1]
    contents = parse_txt_file(input_path)
    text, pattern = contents

    count = pattern_count(text, pattern)
    print("\nNumber of occurrences: {}.\n".format(count))
