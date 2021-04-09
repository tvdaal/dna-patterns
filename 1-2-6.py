"""This module is used for exc. 1.2.6 of the Bioinformatics Specialization.

  Example of how to run:

  python 1-2-6.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Dataset_1.2.6.txt
"""


import sys
import sequence as seq


input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text, pattern = contents

sequence = seq.Sequence(text)
count = sequence.pattern_count(pattern)
print("\nNumber of pattern occurrences:\n\n{}\n".format(count))
