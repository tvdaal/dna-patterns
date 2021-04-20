#!/usr/bin/env python3
"""This script is used to identify DnaA boxes in a given genome.

After providing an input text file in FASTA format containing a genome, the
script identifies G-C skew minima that provide hints as to where the
replication origin might be located. Next, it searches for clumps of patterns
in windows around those minima. Slight mismatches are allowed for, reflecting the
possibility of random nucleotide mutations. Also reverse complements are
included. The patterns that occur most frequently in certain windows are good
candidates for the DnaA boxes of the genome.

  Example of how to run:

  python find_dnaa_box.py <path_to_input_fasta_file>.txt

The various parameters in this script are currently not configurable by the
user. This might change in a future update.
"""


import sequence as seq
import sys
import time


start_time = time.time()
print("\nParsing input file...")
input_path = sys.argv[1]
text = seq.parse_fasta_file(input_path)
sequence = seq.Sequence(text)
print("\nSuccesfully parsed {} into a sequence of length {}.\n".format(input_path, sequence.length))

print("\nDetermining the G-C skew minima...")
skew = sequence.skew_graph()
skew_minima = [str(minimum) for minimum in skew["Skew minima"]]
skew_minima = " ".join(skew_minima)
print("\nLocations of skew minima:\n\n{}\n".format(skew_minima))

print("\nFinding clumps of patterns around the skew minima...")
pattern_length = 9
window_length = 1000
min_occurrences = 6
mismatches = 1
reverse_comp = True
width = int(window_length / 2)
dnaa_boxes = []
for minimum in skew["Skew minima"]:
    min = int(minimum)
    start_pos = (min - width) if (min >= width) else 0
    end_pos = (min + width) if (min <= sequence.length - width) else sequence.length
    window = sequence.sequence[start_pos:end_pos]
    window = seq.Sequence(window)
    clumps = window.find_clumps(
        pattern_length,
        window_length,
        min_occurrences,
        hamming_max=mismatches,
        reverse=reverse_comp,
    )
    dnaa_boxes.extend(clumps["Patterns"])
dnaa_boxes = list(set(dnaa_boxes))
patterns = " ".join(dnaa_boxes)
print("\nA number of {} possible DnaA boxes has been found:\n\n{}\n".format(len(dnaa_boxes), patterns))

elapsed_time = round(time.time() - start_time, 2)
print("\nThis program took", elapsed_time, "seconds to run.\n")
