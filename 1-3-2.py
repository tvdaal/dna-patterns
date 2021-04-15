#!/usr/bin/env python3
"""This module is used for exc. 1.3.2 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-3-2.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-3-2.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]

sequence = seq.Sequence(text)
rc = sequence.reverse_complement()
print("\nReverse complement:\n\n{}\n".format(rc))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
