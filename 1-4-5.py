"""This module is used for exc. 1.4.5 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-4-5.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-4-5.txt
"""


import sys
import time
import sequence as seq


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text, k, L, t = contents

sequence = seq.Sequence(text)
patterns = sequence.find_clumps(int(k), int(L), int(t))
patterns = " ".join(patterns)
print("\nDistinct k-mers that form clumps:\n\n{}\n".format(patterns))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")
