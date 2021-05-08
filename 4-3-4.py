#!/usr/bin/env python3
"""This module is used for exc. 4.3.4 of the Bioinformatics Specialization.

  Example of how to run:

  python 4-3-4.py /Users/tvdaal/Dropbox/Tom/CS/Bioinformatics/Datasets/4-3-4.txt
"""


import sequence as seq
import sys
import time


start_time = time.time()
input_path = sys.argv[1]
contents = seq.parse_txt_file(input_path)
text_1 = contents[0]
text_2 = contents[1]
text_3 = contents[2]
text_lst = contents[3:]

sequences = [seq.Sequence(text) for text in text_lst]
num_iterations = 20
motif_results = seq.repeated_motif_search(
    sequences,
    int(text_1),
    num_iterations,
    gibbs_inner_loop=int(text_3),
)
best_motifs = [motif.sequence for motif in motif_results["Motif matrix"]]
best_motifs = "\n".join(best_motifs)
print("\nThe motif matrix for the collection of {} sequences:\n\n{}\n".format(text_2, best_motifs))
consensus_string = motif_results["Consensus"].sequence
print("The consensus string:\n\n{}\n".format(consensus_string))
print("The motif score:\n\n{}\n".format(motif_results["Score"]))

elapsed_time = round(time.time() - start_time, 2)
print("This program took", elapsed_time, "seconds to run.\n")
