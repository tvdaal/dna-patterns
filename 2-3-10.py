#!/usr/bin/env python3
"""This module is used for exc. 2.3.10 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-3-10.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-3-10.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]

sequence = seq.Sequence(text)
results = sequence.skew_graph()
skew_minima = [str(minimum) for minimum in results["Skew minima"]]
skew_minima = " ".join(skew_minima)
print("\nSkew minima:\n\n{}\n".format(skew_minima))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
