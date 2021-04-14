#!/usr/bin/env python3
"""This module is used for exc. 1.2.6 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-2-6.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-2-6.txt
"""


import sys
import sequence as seq
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text, pattern = contents

sequence = seq.Sequence(text)
count, _ = sequence.pattern_count(pattern)
print("\nNumber of pattern occurrences:\n\n{}\n".format(count))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")
