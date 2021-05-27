#!/usr/bin/env python3
"""This script is used to identify a regulatory motif in a set of sequences.

After providing an input text file containing a list of sequences, the
script runs an algorithm to determine the consensus string. The algorithm is
either a Gibbs search, a randomized search, a greedy search, or a median search.
The greedy and median search algorithms are deterministic and are thus only run
once. The latter is essentially a brute force search for the best solution, so
this algorithm is very slow and can only be used for short motifs. The Gibbs
and randomized search algorithms are fastest and non-deterministic and can be
run repeatedly to find the best result (i.e. the motif matrix with the lowest
motif score). The consensus string represents the regulatory motif.

    Example of how to run:

    python find_motif.py <path_to_input_file>.txt

    For more information on the various optional parameters run:

    python find_motif.py --help
"""


# Standard Python libraries:
from argparse import ArgumentParser, ArgumentTypeError
import time

# Code repository imports:
import sequence as seq


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
        help="Specify the path to the input text file.",
    )
    parser.add_argument(
        "--motif_length",
        type=int,
        default=15,
        help="Specify the length of the motif that will be searched for.",
    )
    parser.add_argument(
        "--algorithm",
        type=str,
        default="gibbs",
        choices=["gibbs", "random", "greedy", "median"],
        help="Specify the algorithm that will be used.",
    )
    parser.add_argument(
        "--num_iterations",
        type=int,
        default=100,
        help="Specify the number of runs for the Gibbs or the randomized search algorithm.",
    )
    parser.add_argument(
        "--gibbs_inner_loop",
        type=int,
        default=1000,
        help="Specify the number of times Gibbs sampling takes place.",
    )
    parser.add_argument(
        "--laplace",
        type=str_to_bool,
        default=True,
        help="Specify whether pseudocounts are added to the profile matrix.",
    )

    return parser


if __name__ == "__main__":
    start_time = time.time()

    print("\nParsing input file...")
    parser = construct_argparser()
    args = parser.parse_args()
    text_lst = seq.parse_txt_file(args.input_path)
    sequences = [seq.Sequence(text) for text in text_lst]
    print("\nSuccesfully parsed {} into a set of {} sequences.\n".format(args.input_path, len(sequences)))

    algorithm = args.algorithm
    print("\nFinding motifs across the sequences with the {} search algorithm...".format(algorithm))
    if algorithm == "median":
        motif_results = seq.median_motif_search(sequences, args.motif_length)
    elif algorithm == "greedy":
        motif_results = seq.greedy_motif_search(
            sequences,
            args.motif_length,
            args.laplace,
        )
    else:
        gibbs = True if (algorithm == "gibbs") else False
        motif_results = seq.repeated_motif_search(
            sequences,
            args.motif_length,
            args.num_iterations,
            gibbs,
            args.gibbs_inner_loop,
            args.laplace,
        )
    best_motifs = [motif.sequence for motif in motif_results["Motif matrix"]]
    best_motifs = "\n".join(best_motifs)
    print("\nThe motif matrix:\n\n{}\n".format(best_motifs))
    consensus_string = motif_results["Consensus"].sequence
    print("The consensus string:\n\n{}\n".format(consensus_string))
    print("The motif score:\n\n{}\n".format(motif_results["Score"]))

    elapsed_time = round(time.time() - start_time, 2)
    print("\nThis program took", elapsed_time, "seconds to run.\n")
