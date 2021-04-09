"""This module is used for exc. 1.4.6 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-4-6.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-4-6.txt
"""


import sys
import time
import sequence as seq


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]
k = 9
L = 500
t = 3

sequence = seq.Sequence(text)
patterns = sequence.find_clumps(k, L, t)
different_patterns = len(patterns)
patterns = " ".join(patterns)
print("\nThere are {} distinct k-mers that form clumps:\n\n{}\n".format(different_patterns, patterns))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")
