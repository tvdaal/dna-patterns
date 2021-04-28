#!/usr/bin/env python3
"""This script is used to identify DnaA boxes in a given genome.

After providing an input text file in FASTA format containing a genome, the
script identifies G-C skew minima that provide hints as to where the
replication origin might be located. Next, it searches for clumps of patterns
in windows around those minima. Slight mismatches are allowed for, reflecting
the possibility of random nucleotide mutations. Also reverse complements are
included. The patterns that occur most frequently in certain windows are good
candidates for the DnaA boxes of the genome.

    Example of how to run:

    python find_dnaa_box.py <path_to_input_fasta_file>.txt

    For more information on the various optional parameters run:

    python find_dnaa_box.py --help
"""


from argparse import ArgumentParser, ArgumentTypeError
import sequence as seq
import sys
import time


def construct_argparser() -> ArgumentParser:
    """Constructs an argument parser for command line use.

    Returns:
        An argument parser.
    """

    def str_to_bool(arg):
        if isinstance(arg, bool):
            return arg
        if arg.lower() in ("yes", "true", "y", "1"):
            return True
        elif arg.lower() in ("no", "false", "n", "0"):
            return False
        else:
            raise ArgumentTypeError("Boolean value expected.")

    parser = ArgumentParser()
    parser.add_argument(
        "input_path",
        type=str,
        help="Specify the path to the input text file in FASTA format.",
    )
    parser.add_argument(
        "--pattern_length",
        type=int,
        default=9,
        help="Specify the path to the input text file in FASTA format.",
    )
    parser.add_argument(
        "--window_length",
        type=int,
        default=1000,
        help="Specify the length of the window in the sequence that is looked at.",
    )
    parser.add_argument(
        "--min_occurrences",
        type=int,
        default=6,
        help="Specify the minimum number of occurrences of the pattern in the given window.",
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=1,
        help="Specify the maximum number of allowed mismatches between similar patterns.",
    )
    parser.add_argument(
        "--reverse",
        type=str_to_bool,
        default=True,
        help="Specify whether reverse complementary patterns should be taken into account as well.",
    )

    return parser


if __name__ == "__main__":
    start_time = time.time()

    print("\nParsing input file...")
    parser = construct_argparser()
    args = parser.parse_args()
    text = seq.parse_fasta_file(args.input_path)
    sequence = seq.Sequence(text)
    print("\nSuccesfully parsed {} into a sequence of length {}.\n".format(args.input_path, sequence.length))

    print("\nDetermining the G-C skew minima...")
    skew = sequence.skew_graph()
    skew_minima = [str(minimum) for minimum in skew["Skew minima"]]
    skew_minima = " ".join(skew_minima)
    print("\nLocations of skew minima:\n\n{}\n".format(skew_minima))

    print("\nFinding clumps of patterns around the skew minima...")
    width = int(args.window_length / 2)
    dnaa_boxes = []
    for minimum in skew["Skew minima"]:
        min = int(minimum)
        start_pos = (min - width) if (min >= width) else 0
        end_pos = (min + width) if (min <= sequence.length - width) else sequence.length
        window = sequence.sequence[start_pos:end_pos]
        window = seq.Sequence(window)
        clumps = window.find_clumps(
            args.pattern_length,
            args.window_length,
            args.min_occurrences,
            hamming_max=args.mismatches,
            reverse=args.reverse,
        )
        dnaa_boxes.extend(clumps["Patterns"])
    dnaa_boxes = list(set(dnaa_boxes))
    patterns = " ".join(dnaa_boxes)
    print("\nA number of {} possible DnaA boxes has been found:\n\n{}\n".format(len(dnaa_boxes), patterns))

    elapsed_time = round(time.time() - start_time, 2)
    print("\nThis program took", elapsed_time, "seconds to run.\n")
