# Introduction
This repository contains tools to analyze DNA sequences. I decided to build it
to support a number of operations and analyses performed in a course on
Bioinformatics, offered by the University of California, San Diego
(see [Coursera](https://www.coursera.org/learn/dna-analysis)). Coding up the
various algorithms made the course much more fun and greatly improved my
understanding of DNA analyses.

# Installation
This repository should work for Python 3.7+. Upon executing
*conda env create -f environment.yml* in the terminal, a virtual conda
environment called 'bio' is created that contains all necessary packages and
dependencies.

# Modules
The main module of this repository is *sequence.py*. It contains functions to
parse an input file with DNA information. Most importantly, it contains the
class Sequence that represents a string of DNA. This class comprises various
methods to perform operations or analyses on DNA sequences. Lastly, the module
contains various algorithms to search for motifs in sets of DNA sequences.

There are two scripts: *find_dnaa_box.py* and *find_motif.py*. Both scripts
make use of the Sequence class. The former script takes a genome as input and
outputs the DnaA box, which is the binding site of the DnaA protein that
initiates DNA replication. The latter script finds regulatory motifs in a
collection of DNA strings. These motifs are short sequences that control the
expression of genes. To identify motifs, one can choose from various algorithms,
some of which are nondeterministic and more suitable for longer motifs.

# Examples
This repository contains two datasets that were provided by the course. The
file *salmonella_enterica.txt* contains the entire genome of the Salmonella
enterica bacterium. The other one, *mtb_dosr.txt*, contains the upstream
regions (with length 250) of 25 genes in Mycobacterium tuberculosis bacterium
(MTB) whose expression is controlled by the dormancy survival regulator (DosR)
transcription factor.

**Example 1:**
Running *python find_dnaa_box.py data/salmonella_enterica.txt* provides
candidates for DnaA boxes in the replication origin of Salmonella enterica
(with the default parameters). Those boxes correspond to the locations where
the DnaA protein binds. This protein activates the initiation of DNA
replication in bacteria.

**Example 2:**
Running *python find_motif.py data/mtb_dosr.txt* provides a motif matrix and
corresponding consensus string for MTB genes that are regulated by DosR (with
the default parameters). The consensus string is a candidate for the motif that
corresponds to the binding site for DosR. This information is useful to
understand when a Tuberculosis infection is active or dormant.

# Improvements
The code can be further improved by implementing proper exception handling to
increase robustness. For the time being, most constraints on input arguments
are implicit (and obvious).
