#!/usr/bin/env python3
"""This module is used for exc. 3.7.1 of the Bioinformatics Specialization.

  Example of how to run:

  python 3-7-1.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/3-7-1.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_lst = contents[1:]

pattern = seq.Sequence(text_1)
sequences = [seq.Sequence(text) for text in text_lst]
total_distance = seq.distance_pattern_strings(sequences, pattern)
print("\nSum of distances between the pattern and each string:\n\n{}\n".format(total_distance))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
