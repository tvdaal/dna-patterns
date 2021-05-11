# bioinformatics
This repository contains tools to analyze DNA sequences. I decided to build it
to support a number of analyses performed in a course on Bioinformatics,
offered by the University of California, San Diego
(see www.coursera.org/learn/dna-analysis).

This repository should work for Python 3.7+. Upon executing
'conda env create -f environment.yml' in the terminal, a virtual conda
environment called 'bio' is created that contains all necessary packages and
dependencies.

The main module of this repository is sequence.py. It contains functions to
parse an input file with DNA information. Most importantly, it contains the
class Sequence that represents a string of DNA. This class comprises various
methods to perform operations or analyses on DNA sequences. Lastly, the module
contains various algorithms to search for motifs in sets of DNA sequences.

There are two scripts: find_dnaa_box.py and find_motif.py. Both scripts make
use of the Sequence class. The former script takes a genome as input and
outputs the DnaA box, which is the binding site of the DnaA protein that
initiates DNA replication. The latter script finds regulatory motifs in a
collection of DNA strings. These motifs are short sequences that control the
expression of genes. To identify motifs, one can choose from various algorithms,
some of which are nondeterministic and more suitable for longer motifs.

The code can be further improved by implementing proper exception handling. For
the time being, most constraints on input arguments are implicit (and obvious).
However, especially for third-party users, robustness of code is essential.
