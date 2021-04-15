#!/usr/bin/env python3
"""This module is used for exc. 1.3.6 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-3-6.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-3-6.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_2 = "CTTGATCAT"

sequence = seq.Sequence(text_1)
pattern = seq.Sequence(text_2)
results = sequence.pattern_count(pattern)
positions_list = results["Positions"]
position_strings = [str(position) for position in positions_list]
positions = " ".join(position_strings)
print("\nStarting positions of the pattern:\n\n{}\n".format(positions))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
