#!/usr/bin/env python3
"""This module is used for exc. 2.7.4 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-7-4.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-7-4.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1, text_2 = contents

sequence = seq.Sequence(text_1)
neighborhood = sequence.neighbors(int(text_2))
neighbors = [str(pattern) for pattern in neighborhood]
neighbors = " ".join(neighbors)
print("\nNeighborhood:\n\n{}\n".format(neighbors))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
