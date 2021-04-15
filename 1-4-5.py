#!/usr/bin/env python3
"""This module is used for exc. 1.4.5 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-4-5.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-4-5.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1, text_2, text_3, text_4 = contents

sequence = seq.Sequence(text_1)
results = sequence.find_clumps(int(text_2), int(text_3), int(text_4))
patterns = " ".join(results["Patterns"])
print("\nDistinct patterns that form clumps in the given sequence:\n\n{}\n".format(patterns))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
