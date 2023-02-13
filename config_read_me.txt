# Read me/config file for MB_PRIME.
# Author: Ben Tudor Price
# Email: benjamin.tudorprice@mail.utoronto.ca
# Version History
#   08/04/2022: V1.0
#       - Initial release.

# This program was initially designed for two step PCR, which requres the
# adapter region specified below. However, one can still add the entirety of the
# 5' region of the primer in lieu of the adapter, which will result in valid
# one step PCR primers. This program uses rudementary base counts as its core
# metric when evaluating dimer formation potential, so the use of external
# dimer evaluation metrics is reccomended.
#

# === Primer Structure Overview ===
# 5'[Adapter Sequence]-[Heterogeneity Spacer]-[Binding Sequence]- 3'

# This program works to design heterogeneity spacers that do not introduce
# additional dimer structures into the primers. Here, we'll need the sequences
# to which you'd like your primers to bind, 5' - 3'.

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

*SHOW_SPACER_MENU = True

# This program is capable of multithreading to decrease runtime. Enter the
# number of cores you'd like to use below.

*NUM_CORES = 1

# The program is also capable of increasing the number of heterogeneity spacers
# sampled, and thus the probability of finding better primer. If this is really
# important to you, increasing the rigour will marginally increase primer quality
# at the cost of an exponential runtime increase.

*RIGOUR = 1


# The number of sets to generate. Does not affect runtime.
*NUM_SETS_TO_GEN = 3












