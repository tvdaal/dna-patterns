#!/usr/bin/env python3
"""This module is used for exc. 2.6.2 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-6-2.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-6-2.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = "".join(contents)

sequence = seq.Sequence(text)
results = sequence.skew_graph()
skew_minima = [str(minimum) for minimum in results["Skew minima"]]
skew_minima = " ".join(skew_minima)
print("\nSkew minima:\n\n{}\n".format(skew_minima))

patterns = sequence.frequent_words(9, hamming_max=1)
patterns = " ".join(patterns)
print("\nMost frequent (approximate) patterns:\n\n{}\n".format(patterns))

results = sequence.find_clumps(9, 500, 4, hamming_max=1)
patterns = " ".join(results["Patterns"])
print("\nDistinct patterns that form clumps in the given sequence:\n\n{}\n".format(patterns))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
