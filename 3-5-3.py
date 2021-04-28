#!/usr/bin/env python3
"""This module is used for exc. 3.5.3 of the Bioinformatics Specialization.

  Example of how to run:

  python 3-5-3.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/3-5-3.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_2 = contents[1]
pattern_length = int(text_2)
text_3 = contents[2:2+pattern_length]
text_4 = contents[2+pattern_length:2+2*pattern_length]
text_5 = contents[2+2*pattern_length:2+3*pattern_length]
text_6 = contents[2+3*pattern_length:2+4*pattern_length]

profile_matrix = {
    "A": [float(score) for score in text_3],
    "C": [float(score) for score in text_4],
    "G": [float(score) for score in text_5],
    "T": [float(score) for score in text_6],
}
sequence = seq.Sequence(text_1)
most_prob_string = sequence.most_probable_string(pattern_length, profile_matrix)
print("\nProfile-most probable string:\n\n{}\n".format(most_prob_string))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
