# Read me/config file for MB_PRIME.
#
# The three executables included in this folder are intended to help at each
# stage of the primer design process for metabarcoding. Each of these has a
# detailed description of its purpose and usage instructions provided below.
#
# Written for python 3.10
#
# Author: Ben Tudor Price
# Email: benjamin.tudorprice@mail.utoronto.ca
# Version History
#   08/04/2022: V1.0
#       - Initial beta version.

# This file contains a list of parameters that will be used at each stage in the
# program. All of the variables under any one of the program headers must have
# some value in order to run. Values that must be stated are marked with a *. Any
# line that begins with a "#" will be ignored during parsing, so you are free
# to delete and add them as you wish.
#
# Ensure that you have python3.10 installed, else you may get crypic errors
# involving match blocks.
#
# To install the required packages, run the following on linux.
# > pip install pipreqs
# > pip install -r requirements.txt


# === Primer Structure Overview ===
# 5'[Adapter Sequence]-[Indexing Sequence]-[Heterogeneity Spacer]-[Binding Sequence]- 3'


# === Program 1: analyse_alignment.py ===

*MSA_FILEPATH = C:\Users\bfern\PycharmProjects\mb_prime\test_files\test_objects\alignments\alex_test_align.fas

*GRAPH_OUT_FILEPATH = DIR

*STORE_IMAGE = False
*DISPLAY_GRAPH = True

# === Program 2: scan_for_binding_regions.py ===


*FORWARD_ADAPTERS =
@BEGIN_FASTA
> 1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
@END_FASTA

*REVERSE_ADAPTERS =
@BEGIN_FASTA
>1
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
@END_FASTA

*FORWARD_BINDING_START = 100
*FORWARD_BINDING_END = 300

*REVERSE_BINDING_START = 400
*REVERSE_BINDING_END = 600

*FORWARD_BINDING_LENGTHS = 15
*REVERSE_BINDING_LENGTHS = 15

*AMPLICON_LENGTH_MIN = 100
*AMPLICON_LENGTH_MAX = 200

*CONSERVATION_WEIGHT = 6
*DIMER_WEIGHT = 1


=== Program 3: heterogeneity_spacer_gen.py ===

*FORWARD_BINDING_SEQ = #str#
*REVERSE_BINDING_SEQ = #str#

*FORWARD_ADAPTER = #str#
*REVERSE_ADAPTER = #str#

*HETEROGENEITY_REGION_LENGTH = #int#

*OUTPUT_FASTA = DIR

*SHOW_SPACER_MENU = False

*NUM_CORES = 1

*RIGOUR = 3













