#!/usr/bin/env python3
"""This module is used for exc. 2.4.3 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-4-3.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-4-3.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1, text_2 = contents

sequence_1 = seq.Sequence(text_1)
sequence_2 = seq.Sequence(text_2)
hamming_distance = sequence_1.hamming_distance(sequence_2)
print("\nThe Hamming distance between the two sequences:\n\n{}\n".format(hamming_distance))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
