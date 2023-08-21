[![DOI](https://zenodo.org/badge/601367802.svg)](https://zenodo.org/badge/latestdoi/601367802)
# hetero_spacer_gen

hetero_spacer_gen is a program than enables the design of heterogeneity spacers that introduce nucleotide diversity into low diversity libraries.

# Installation

This program requires python3.10. To install the required packages, run the following.
> pip install pipreqs
> 
> pip install -r requirements.txt


# Usage
Before executing the program, enter the other the neccesary primer components in the 'config_read_me.txt' file. In order to start the program, execute:

>python3 run.py

This will run the program using the configuration parameters in the 'config_read_me.txt' file stored in the project directory and output primer sets in a .fasta in the specified location.

# Example Use Case

We'd like to perform one-step PCR to amplify a region of DNA, which we will then sequence using an Illumina platform. The final structure such a primer is as follows:

5' - Grafting - Index - Sequencing Primer - Heterogeneity Spacer - Binding Region - 3'

Where the binding region is responsible for binding the target DNA molecule. Note that in the config file we refer to the "Grafting - Index - Sequencing Primer" portion of the primer as an adapter, which is a quirk of this program being initially designed for two-step PCR.

Now, we need to copy our primer components into the config file (5'-3').

```
*FORWARD_BINDING_SEQ = ATCGAGAGCTAGCTAGTCGATCGGC
*REVERSE_BINDING_SEQ = GGATCGGATCGGGATGGGCTATACC

*FORWARD_ADAPTER = ATCGATCGATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGAGTGTGTTCGGGAGGGATC
*REVERSE_ADAPTER = ATCGATGATGTACGTGACTAGGATGCTGTAGCTAGTCGATCGATGTTTATGCTAGAGCGATGTGCTGATGAA

*HETEROGENEITY_REGION_LENGTH = 12
```
The heterogeneity length specifies the number of bases that should be perfectly diverse - this can often be accomplished using spacers that are shorter than this length.

Now we run:

> python3.10 run.py

We can then select the pairing of heterogeneity spacer length (option 1 is usually the best). The program then performs some small optimisation steps which will slightly reduce dimerisation potential for the given primers. Our primer set may look similar to the one below, with the heterogeneity spacers in lower case.
```
>Forward #1
ATCGATCGATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGAGTGTGTTCGGGAGGGATCATCGAGAGCTAGCTAGTCGATCGGC
>Reverse #1
ATCGATGATGTACGTGACTAGGATGCTGTAGCTAGTCGATCGATGTTTATGCTAGAGCGATGTGCTGATGAAgcGGATCGGATCGGGATGGGCTATACC
>Forward #2
ATCGATCGATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGAGTGTGTTCGGGAGGGATCtgtATCGAGAGCTAGCTAGTCGATCGGC
>Reverse #2
ATCGATGATGTACGTGACTAGGATGCTGTAGCTAGTCGATCGATGTTTATGCTAGAGCGATGTGCTGATGAAatacGGATCGGATCGGGATGGGCTATACC
>Forward #3
ATCGATCGATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGAGTGTGTTCGGGAGGGATCgagcgtctagcATCGAGAGCTAGCTAGTCGATCGGC
>Reverse #3
ATCGATGATGTACGTGACTAGGATGCTGTAGCTAGTCGATCGATGTTTATGCTAGAGCGATGTGCTGATGAAtacttctaatcGGATCGGATCGGGATGGGCTATACC
>Forward #4
ATCGATCGATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGAGTGTGTTCGGGAGGGATCccatcatctcttATCGAGAGCTAGCTAGTCGATCGGC
>Reverse #4
ATCGATGATGTACGTGACTAGGATGCTGTAGCTAGTCGATCGATGTTTATGCTAGAGCGATGTGCTGATGAAcgtacagctcatGGATCGGATCGGGATGGGCTATACC
```
