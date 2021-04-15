#!/usr/bin/env python3
"""This module is used for exc. 1.4.6 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-4-6.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-4-6.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_2 = 9
text_3 = 500
text_4 = 3

sequence = seq.Sequence(text_1)
results = sequence.find_clumps(text_2, text_3, text_4)
different_patterns = results["Number"]
patterns = " ".join(results["Patterns"])
print("\nWe have found the following {} distinct patterns that form clumps:\n\n{}\n".format(different_patterns, patterns))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
