"""This module is used for exc. 1.2.13 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-2-13.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Dataset_1.2.13.txt
"""


import sys
import sequence as seq


input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text, k = contents

sequence = seq.Sequence(text)
patterns = sequence.frequent_words(int(k))
print("\nMost frequent patterns:\n\n{}\n".format(patterns))
