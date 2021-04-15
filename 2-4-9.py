#!/usr/bin/env python3
"""This module is used for exc. 2.4.9 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-4-9.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-4-9.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1, text_2, text_3 = contents

sequence = seq.Sequence(text_1)
patterns = sequence.frequent_words(int(text_2), hamming_max=int(text_3))
patterns = " ".join(patterns)
print("\nMost frequent (approximate) patterns:\n\n{}\n".format(patterns))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
