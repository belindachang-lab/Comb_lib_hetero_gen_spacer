# Read me/config file for MB_PRIME.
#
# The three executables included in this folder are intended to help at each
# stage of the primer design process for metabarcoding. Each of these has a
# detailed description of its purpose and usage instructions provided below.
#
#
#
# Author: Ben Tudor Price
# Email: benjamin.tudorprice@mail.utoronto.ca
# Version History
#   08/04/2022: V1.0
#       - Initial release.

# This file contains a list of parameters that will be used at each stage in the
# program. All of the variables under any one of the program headers must have
# some value in order to run. Values that must be stated are marked with a *. Any
# line that begins with a "#" will be ignored during parsing, so you are free
# to delete and add them as you wish.


# === Primer Structure Overview ===
# 5'[Adapter Sequence]-[Indexing Sequence]-[Heterogeneity Spacer]-[Binding Sequence]- 3'


# === Program 1: analyse_alignment.py ===

# This program is meant to give the user a general understanding of the
# alignment that they'd like to design metabarcoding primers for. Given an MSA,
# this program will display a graph containing several properties useful when
# deciding where to place our primers.
#
# Conservation:
#     This is important, as primers targeting less conserved regions will not
#     bind all sequences in the alignment with equal affinity, potentially
#     resulting in a failure to amplify those sequences.
#
# Sequences Missed:
#     Regions containing indels are hard to target with primers, as misalignment
#     between the primer and some subset of the sequences greatly reduces
#     primer affinity for the target sequences. This metric shows what % of
#     sequences will be lost due to this, at a minimum.
#

# The first parameter required is the filepath of the MSA to be used. At the,
# moment, only '.fsa' formatted sequences are accepted.

*MSA_FILEPATH = #str#

# You can also specify an outpath for an image of the graph, and whether you
# want to display/store the graph. HEAD just means that the repo path will be
# used.

*GRAPH_OUT_FILEPATH = DIR
*WINDOW_SIZE = 10
*STORE_IMAGE = False
*DISPLAY_GRAPH = True

# === Program 2: scan_for_binding_regions.py ===

# Now that we have a rough idea as to where we'd like our primers to bind, we
# can optimise the binding regions we select. Currently, optimisation includes
# the following criteria.
#
# Conservation:
#       Binding regions with the highest average conservation will be preferred.
#
# Dimerisation Potential:
#       The program will avoid selecting regions that have high complementarity
#       to the adapter sequences, and themselves.

# First, we'll need the sequences of the adapter sequences you plan to use in
# the regions below in a fasta format. Ensure there are an equal number of each.
# If you have indexing sequences that you'd like included, concatenate them
# to these adapters.

# Example:
# *FORWARD_ADAPTERS @
# @BEGIN_FASTA
# >f1
# ATGCTAGCATGCATGTGATCGTAGCGCCGCGG
# >f2
# ATGCTAGCATGCATGTGATCGTAGCGCCGCGG
# @END_FASTA

*FORWARD_ADAPTERS =
@BEGIN_FASTA
@END_FASTA

*REVERSE_ADAPTERS =
@BEGIN_FASTA
@END_FASTA

# Below, specify the regions you would like to search for primer binding
# sites (inclusive). These should be the position to which the most 5' base in
# each primer will bind within the consensus sequence of your alignment.

*FORWARD_BINDING_START = #int#
*FORWARD_BINDING_END = #int#

*REVERSE_BINDING_START = #int#
*REVERSE_BINDING_END = #int#

# Next we need the allowable lengths of the binding sequence (see primer
# diagram).

# EXAMPLE
#FORWARD_BINDING_LENGTHS = 12, 13, 14

*FORWARD_BINDING_LENGTHS = #int, int, ...#
*REVERSE_BINDING_LENGTHS = #int, int, ...#

# Finally, we need the allowable amplicon lengths. These should be compatible
# with the binding regions you've specified (e.g. if your only forward site is at
# base 300, and your only reverse site is 400, then you must have an amplicon of
# at least size 101).

*AMPLICON_LENGTH_MIN = #int#
*AMPLICON_LENGTH_MAX = #int#

# Optionally, you can specify the weighting of the criteria. Turn to zero to
# ignore completely.

*CONSERVATION_WEIGHT = 1
*DIMER_WEIGHT = 2


=== Program 3: heterogeneity_spacer_gen.py ===

# This program works to design heterogeneity spacers that do not introduce
# additional dimer structures into the primers. Here, we'll need the binding
# sequences that you obtained from the previous program. 5'-3'

*FORWARD_BINDING_SEQ = #str#
*REVERSE_BINDING_SEQ = #str#

# You'll also need to specify the adapters you'd like to use. 5'-3'
*FORWARD_ADAPTER = #str#
*REVERSE_ADAPTER = #str#

# Below, enter the number of bases you'd like to ensure heterogeneity for.

*HETEROGENEITY_REGION_LENGTH = #int#

# Finally, enter the location of the output fasta location.

*OUTPUT_FASTA = DIR

# By default, the program will select the shortest set of heterogeneity spacers
# possible, you can prompt the program to display a list of these if you'd like.

*SHOW_SPACER_MENU = False

# This program is capable of multithreading to decrease runtime. Enter the
# number of cores you'd like to use below.

*NUM_CORES = 1

# The program is also capable of increasing the number of heterogeneity spacers
# sampled, and thus the probability of finding better primer. If this is really
# important to you, increasing the rigour will marginally increase primers at
# the cost of an exponential runtime increase.

*RIGOUR = 3

# The number of sets to generate. Does not affect runtime.
*NUM_SETS_TO_GEN = 3













