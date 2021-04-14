#!/usr/bin/env python3
"""This module is used for exc. 2.3.8 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-3-8.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-3-8.txt
"""


import sys
import sequence as seq
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]

sequence = seq.Sequence(text)
skew_scores, _ = sequence.skew_graph()
skew_scores = [str(score) for score in skew_scores]
skew_scores = " ".join(skew_scores)
print("\nSkew scores:\n\n{}\n".format(skew_scores))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")
