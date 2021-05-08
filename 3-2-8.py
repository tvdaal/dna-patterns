#!/usr/bin/env python3
"""This module is used for exc. 3.2.8 of the Bioinformatics Specialization.

  Example of how to run:

  python 3-2-8.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/3-2-8.txt
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

sequences = [seq.Sequence(text) for text in text_lst]
patterns = seq.find_motifs(sequences, int(text_1), int(text_2))
motifs = [str(pattern) for pattern in patterns]
num_motifs = len(motifs)
motifs = " ".join(motifs)
print("\nA number of {} of regulatory motifs were found:\n\n{}\n".format(num_motifs, motifs))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
