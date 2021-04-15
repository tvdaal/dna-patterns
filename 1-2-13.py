#!/usr/bin/env python3
"""This module is used for exc. 1.2.13 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-2-13.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-2-13.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1, text_2 = contents

sequence = seq.Sequence(text_1)
patterns = sequence.frequent_words(int(text_2))
patterns = " ".join(patterns)
print("\nMost frequent patterns:\n\n{}\n".format(patterns))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
