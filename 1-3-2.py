"""This module is used for exc. 1.3.2 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-3-2.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Dataset_1.3.2.txt
"""


import sys
import sequence as seq


input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]

sequence = seq.Sequence(text)
rc = sequence.reverse_complement()
print("\nReverse complement:\n\n{}\n".format(rc))
