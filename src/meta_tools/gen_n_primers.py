import os

from src.hetero_spacer_generator.spacer_generator.spacer_filters import \
    SortForPairwise
from time import time
from src.hetero_spacer_generator.primer_tools import MBPrimerBuilder, MaxInt, \
    PrimerSet
from typing import TypeVar
import random
from src.hetero_spacer_generator.spacer_generator.hetero_spacer_generator import \
    HeteroGen
from pathlib import Path
from src.meta_tools.template_sequences import *
from src.hetero_spacer_generator.defaults import V

DESKTOP = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop')

# Defaults for generation of primers.
RIGOUR = 1
# The number of sets generated for binding sequence - adapter pair.
NUM_SETS = 1
# Output file for primers.
OUT_FILE = DESKTOP / "Demo\\Output Primers"
# Print additional information
VERBOSE = True
# Format the output in fasta format. Required for TF analysis.
DO_FASTA = True
# Starting adapters and binding sequences. See available in
# meta_tools.template_sequences.py
GEN_SET = STANDARD_SET
# Whether output files should actually be made.
GEN_OUTPUT = True
# Whether you'd like seperate files for each binding seq - adapter pair.
MULTI_FILE = False
# Whether you'd like to have a csv file containing the generated scores returned
# A CSV will be created for each call to output_primer_sets
CSV = True
# The number of processes to use
NUM_PROCS = os.cpu_count()
# Whether the heterogeneity spacings used should be variable.
VARY_SPACERS = False

# Size of heterogeneity region
NUM_HETERO = 12
# Maximum spacer length
MAX_SPACER_LENGTH = 12


# I've left the above variables as default values in output_primer_sets, in case
# you'd like to vary them, as shown below.

# Put whatever you'd like to do in here.
def main():
    rand = MaxInt(False)
    rand_rep = 1000
    high = 5
    high_rep = 1
    gen_set = MITOFISH

    # Don't need multiprocessing here.
    global NUM_PROCS
    NUM_PROCS = 1

    filep = DESKTOP / 'Alex Sets' / (gen_set.name + " Random")
    safe_make(filep)
    run_name = gen_set.name
    #output_primer_sets(gen_set, n=rand_rep, rigour=rand,
    #                       filepath=filep, filename=run_name)
    NUM_PROCS = os.cpu_count()

    filep = DESKTOP / 'Alex Sets' / (gen_set.name + " High Rigour 3")
    safe_make(filep)
    run_name = gen_set.name
    output_primer_sets(gen_set, n=high_rep, rigour=-3,
                       filepath=filep, filename=run_name)


def safe_make(path: Path) -> None:
    """Makes the given <path> directory if it does not already exist."""
    if not path.exists():
        path.mkdir()
    else:
        if VERBOSE:
            print(str(path), 'already exists.')


class PrimerGen:
    """Generates sets of primers."""
    _set_list: List[PrimerSet]
    _fresh_sets: bool
    _csv_str: str

    def __init__(self):
        """Initialises an empty primer gen class."""
        self._set_list = []
        self._fresh_sets = False
        self._csv_str = ''

    def get_sets(self) -> List[PrimerSet]:
        """Returns a *newly generated* list of primer sets. Raises an error if
        the sets have been fetched before."""
        if not self._fresh_sets:
            raise AttributeError("Same set has been retreived more than once.")
        return self._set_list

    def get_csv_header(self) -> str:
        """Returns a csv header compatible with any CSVs output by this
        funciton."""
        return 'rigour, ' + SortForPairwise.get_csv_formatting_str() + '\n'

    def get_csv_strs(self) -> str:
        """Returns a string representation of the csv strings accumulated since
        the last call to this method."""
        rtrn = self._csv_str
        self._csv_str = ''
        return rtrn

    def gen_n_primer_sets(self, for_adapter: Seq, rev_adapter: Seq,
                          for_binding: Seq, rev_binding: Seq, n: int,
                          rigour: int, verbose: bool = V) -> None:
        """Returns a list of <n> sets of primers, generated using the given
        components, with <rigour>. Seperates each set with <sep>. Will print
        additional information to console iff <V>. Will format data in fasta format
        iff fasta. Will shift the index of the printed primer by <ind> * <n> bases.
        === Parameters ===
        for_adapter: 5' - 3'
            The forward adapter sequenced used to generate the primer sets.
        rev_adapter: 5' - 3'
            The reverse adapter sequenced used to generate the primer sets.
        for_binding: 5' - 3'
            The forward binding sequenced used to generate the primer sets.
        rev_binding: 5' - 3'
            The reverse binding sequenced used to generate the primer sets.
        n:
            The number of primer sets to generate. Each set uses a different spacer
            combo by default.
        rigour:
            The rigour with which to generate the files.
        """

        hg = HeteroGen(rigour=rigour,
                       num_hetero=NUM_HETERO,
                       max_spacer_length=MAX_SPACER_LENGTH)
        rg = hg.get_primer_gen()
        sf = rg.get_spacer_sorter()

        self._set_list = []

        incomplete_forward_primer = MBPrimerBuilder(adapter_seq=for_adapter,
                                                    binding_seq=for_binding)
        incomplete_reverse_primer = MBPrimerBuilder(adapter_seq=rev_adapter,
                                                    binding_seq=rev_binding)

        # Choose n random spacers from the ten best.
        for_spacers = hg.get_all_spacer_combos(
            incomplete_forward_primer._binding_seq)
        rev_spacers = hg.get_all_spacer_combos(
            incomplete_reverse_primer._binding_seq)

        # Select best spacers
        for_spacer = for_spacers[0]
        rev_spacer = rev_spacers[0]

        for set_num in range(n):

            if VARY_SPACERS:
                for_spacer = random.choice(for_spacers[0:4])
                rev_spacer = random.choice(rev_spacers[0:4])

            number_to_return = 1

            if verbose:
                print(for_spacer)
                print(rev_spacer)
                print("Generating primer set", str(set_num + 1), 'of', str(n),
                      "@rigour", str(rigour) + '.')

            # Generate the sets.
            t0 = time()
            psets = rg.get_hetero_seqs(incomplete_forward_primer,
                                       incomplete_reverse_primer,
                                       for_spacer, rev_spacer, number_to_return,
                                       NUM_PROCS)
            runtime = time() - t0

            if verbose:
                print("Completed set #{num:d} at rigour {rig:d} in {time:.2f}"
                      " seconds.".format(time=runtime, num=set_num + 1,
                                         rig=rigour))

            self._set_list.extend(psets)
            self._fresh_sets = True
            # Append to growing csv string. Convert tup to str. Add rigour to str.
            self._csv_str += str(rigour) + ', ' + str(sf.get_csv())[1:-1] + '\n'


def primer_sets_to_str(primer_sets: List[PrimerSet], fasta: bool, ind: int,
                       sep: str) -> str:
    """Coverts the given PrimerSets into a string formatted as a fasta iff fasta
    else plaintext. Uses <ind> * len(primer_sets) as the starting index for the
    primer names. Divides each set of primers by <sep>."""
    primer_str = ''
    for i, p_set in enumerate(primer_sets):
        if fasta:
            primer_str += p_set.get_fasta_seqs(i + ind * len(primer_sets))
        else:
            primer_str += p_set.get_plain_seqs()
        primer_str += sep

    return primer_str


T = TypeVar('T')


def make_iter(obj: T) -> Union[T, Tuple[T]]:
    """Makes the object iterable if it is not already."""
    try:
        iter(obj)
        return obj
    except TypeError:
        return [obj]


def output_primer_sets(gen_set: GenSet, n: int,
                       rigour: Union[int, List[int]] = RIGOUR,
                       filepath: Path = OUT_FILE,
                       sep: str = '\n', verbose: bool = VERBOSE,
                       fasta: bool = DO_FASTA,
                       multi_file=MULTI_FILE,
                       filename: str = "Primer Set") -> None:
    """Runs gen_n_primers for each of the given binding sequences with the given
     parameters. Stores the output in text files in <filepath>. Will place
     all sets into a single file iff multi_file, othwise will create output file
     for each pair of forward and reverse primers."""
    p_gen = PrimerGen()

    if verbose:
        print("Generating sets for", filename, "along rigours:",
              str(rigour)[1:-1])

    rigour = make_iter(rigour)
    if fasta:
        ft = '.fst'
    else:
        ft = '.txt'

    for rig in rigour:
        for_bindings, rev_bindings, for_adapters, rev_adapters = gen_set.unpack()
        out_str = ''
        ind = 0
        for i in range(min(len(for_bindings), len(rev_bindings))):
            for_binding = for_bindings[i]
            rev_binding = rev_bindings[i]
            for_adapter = rev_adapters[i]
            rev_adapter = for_adapters[i]

            if not multi_file:
                ind = i

            # Get some primer sets as a string.
            p_gen.gen_n_primer_sets(for_adapter, rev_adapter, for_binding,
                                    rev_binding, n, rig)
            sets = p_gen.get_sets()
            s = primer_sets_to_str(sets, fasta, ind, sep)
            txt_path = filepath / "{filename}@set{num:d}@rig{r}{ft}".format(
                filename=filename, num=i + 1, r=rig, ft=ft)
            # Store string in file.
            if multi_file and GEN_OUTPUT:
                with open(txt_path, "w") as my_file:
                    my_file.write(s)
            else:
                out_str += s

        if not multi_file and GEN_OUTPUT:
            txt_path = filepath / "{filename}@rig{r}{ft}".format(
                filename=filename, r=rig, ft=ft)

            if verbose:
                print(filename + " complete")
            with open(txt_path, "w") as my_file:
                my_file.write(out_str)

    if CSV:
        csv_path = filepath / (filename + ' scores.csv')
        with open(csv_path, "w") as my_file:
            my_file.write(p_gen.get_csv_header())
            my_file.write(p_gen.get_csv_strs())


if __name__ == '__main__':
    main()
