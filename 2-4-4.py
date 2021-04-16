#!/usr/bin/env python3
"""This module is used for exc. 2.4.4 of the Bioinformatics Specialization.

  Example of how to run:

  python 2-4-4.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/2-4-4.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_2, text_1, text_3 = contents

sequence = seq.Sequence(text_1)
pattern = seq.Sequence(text_2)
results = sequence.pattern_count(pattern, hamming_max=int(text_3))
positions_list = results["Positions"]
position_strings = [str(position) for position in positions_list]
positions = " ".join(position_strings)
print("\nNumber of pattern occurrences:\n\n{}".format(results["Count"]))
print("\nStarting positions of the given pattern:\n\n{}\n".format(positions))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
