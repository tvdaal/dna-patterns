"""This module is used for exc. 1.3.5 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-3-5.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/1-3-5.txt
"""


import sys
import time
import sequence as seq


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
pattern, text = contents

sequence = seq.Sequence(text)
_, positions_list = sequence.pattern_count(pattern)
position_strings = [str(position) for position in positions_list]
positions = " ".join(position_strings)
print("\nStarting positions of the pattern:\n\n{}\n".format(positions))

elapsed_time = round(time.time() - start_time)
print("This program took", elapsed_time, "seconds to run.\n")