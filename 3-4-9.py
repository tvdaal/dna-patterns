#!/usr/bin/env python3
"""This module is used for exc. 3.4.9 of the Bioinformatics Specialization.

  Example of how to run:

  python 3-4-9.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/3-4-9.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text = contents[0]
text_lst = contents[1:]

sequences = [seq.Sequence(text) for text in text_lst]
motif_results = seq.median_motif_search(sequences, int(text))
best_motifs = [motif.sequence for motif in motif_results["Motif matrix"]]
best_motifs = "\n".join(best_motifs)
print("\nThe motif matrix for the collection of {} sequences:\n\n{}\n".format(len(text_lst), best_motifs))
consensus_string = motif_results["Consensus"].sequence
print("The median (consensus) string:\n\n{}\n".format(consensus_string))
print("The motif score:\n\n{}\n".format(motif_results["Score"]))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
