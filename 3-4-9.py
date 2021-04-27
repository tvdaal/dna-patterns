#!/usr/bin/env python3
"""This module is used for exc. 3.4.9 of the Bioinformatics Specialization.

  Example of how to run:

  python 3-4-9.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/3-4-9.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_2 = contents[1]
text_lst = contents[2:]

sequence = seq.Sequence(text_2)
sequences = [seq.Sequence(text) for text in text_lst]
median_string = sequence.median_string(sequences, int(text_1))
print("\nMedian string:\n\n{}\n".format(median_string))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
