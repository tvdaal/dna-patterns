#!/usr/bin/env python3
"""This module is used for exc. 2.3.10 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-3-10.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-3-10.txt
"""


import sys
import sequence as seq
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]

sequence = seq.Sequence(text)
_, skew_minima = sequence.skew_graph()
print(skew_minima)
skew_minima = [str(minimum) for minimum in skew_minima]
skew_minima = " ".join(skew_minima)
print("\nSkew minima:\n\n{}\n".format(skew_minima))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")
