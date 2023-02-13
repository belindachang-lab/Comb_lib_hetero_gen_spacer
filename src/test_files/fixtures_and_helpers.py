from pathlib import Path
from typing import List, Tuple, TypeVar
import random
from src.hetero_spacer_generator.spacer_generator import SpacerAlignmentGen
from src.hetero_spacer_generator import MBPrimerBuilder
from src.hetero_spacer_generator.spacer_generator.random_spacer_generator import \
    RandomSpacerGen
from src import hetero_spacer_generator as rsg
from Bio.Seq import Seq
from src.hetero_spacer_generator.spacer_generator.spacer_filters import SpacerSet
import os

from seq_alignment_analyser.align import MSA

TEST_PATH = Path(os.path.dirname(__file__))
TEST_OBJECTS_PATH = TEST_PATH / 'test_objects'
ALIGNMENTS_PATH = TEST_OBJECTS_PATH / 'alignments'
TEST_OUTPUT_PATH = TEST_PATH / 'test_output'
READ_MES_PATH = TEST_OBJECTS_PATH / 'read_mes'

LARGE_SAMPLE_SIZE = 100
MED_SAMPLE_SIZE = 50
SMALL_SAMPLE_SIZE = 10

SMALL_SEQ_LEN = 2


def get_msa(filename: str) -> MSA:
    """Returns an MSA parsed from the given <filename> in the alignments
    file."""
    return MSA(ALIGNMENTS_PATH / (filename + '.fas'))

def ensure_hetero_bases(bases: List[str]) -> bool:
    for i in range(len(bases)):
        for j in range(i + 1, len(bases)):
            if bases[i] == bases[j]:
                return False
    return True


def ensure_hetero_seq_arr(seq_arr: [List[List[str]]], num_hetero: int) -> bool:
    """Returns true iff the sequences contained in <seq_arr> are
    heterogeneous for the first <num_hetero> bases."""
    for column in range(min(num_hetero, len(seq_arr[0]) - 1)):
        bases = []
        for row in range(4):
            if seq_arr[row][column].isalpha():
                bases.append(seq_arr[row][column])
        if not ensure_hetero_bases(bases):
            return False
    return True


def ensure_valid_spacers(incomplete_primer: MBPrimerBuilder or Seq,
                         spacers: Tuple[int, int, int, int],
                         num_hetero) -> bool:
    """Ensures that the spacers produced by the given spacers will be valid. If
    """
    if type(incomplete_primer) == MBPrimerBuilder:
        seq_arr = rsg.gen_sequence_array(incomplete_primer.get_binding_seq(),
                                         spacers)
    else:
        seq_arr = rsg.gen_sequence_array(incomplete_primer,
                                         spacers)
    return ensure_hetero_seq_arr(seq_arr, num_hetero)

def gen_random_seq(length: int, GC: float = 0.5) -> Seq:
    """Generates a random Seq of <length>"""
    AT = 1 - GC
    bases = ['A', 'T', 'C', 'G']

    return Seq(''.join(random.choices(bases,
                                      weights=[AT / 2, AT / 2, GC / 2, GC / 2],
                              k=length)))


def gen_random_seq_degen(length: int) -> Seq:
    """Generates a random Seq of <length>"""
    bases = ['A', 'T', 'C', 'G', 'W', 'S',
             'M', 'K', 'R', 'Y', 'B', 'D',
             'H', 'V', 'N']
    seq_str = ''
    for _ in range(length):
        seq_str += random.choice(bases)

    return Seq(seq_str)


def gen_random_spacers(incomplete_forward_primer: MBPrimerBuilder,
                       incomplete_reverse_primer: MBPrimerBuilder,
                       forward_spacer_seqs: List[SpacerSet],
                       reverse_spacer_seqs: List[SpacerSet],
                       forward_spacer: Tuple[int, int, int, int],
                       reverse_spacer: Tuple[int, int, int, int],
                       num_to_generate: int) -> None:
    """Empties <forward_spacer_seqs> and <reverse_spacer_seqs>, refilling them
    with <num_to_generate> randomly generated sets of heterogeneity spacers.
     Generated spacers with be based upon <incomplete_forward_primer> and
     <incomplete_reverse_primer>."""
    forward_spacer_seqs.clear()
    reverse_spacer_seqs.clear()
    forward_spacers = rsg.gen_hetero_set(incomplete_forward_primer,
                                         forward_spacer, num_to_generate)
    forward_spacer_seqs.extend(forward_spacers)
    reverse_spacers = rsg.gen_hetero_set(incomplete_reverse_primer,
                                         reverse_spacer, num_to_generate)
    reverse_spacer_seqs.extend(reverse_spacers)


ADAPTER_LEN = 12
INDEX_LEN = 4
BINDING_MAX = 30
BINDING_MIN = 12


def gen_incomplete_primer(binding_len: int = 12):
    """Randomly generates an incomplete primer with a binding seq with len
    <binding_len> and blank heterogeneity regions. Sets the adapter and
    indexing regions to random sequences of default len."""
    primer = MBPrimerBuilder()
    primer.set_adapter_seq(gen_random_seq(ADAPTER_LEN))
    primer.set_binding_seq(gen_random_seq(binding_len))
    primer.set_index_seq(gen_random_seq(INDEX_LEN))
    return primer


def gen_rand_spacer_aligns(num_primers: int, max_val: int) \
        -> List[Tuple[int, int, int, int]]:
    """Generates a set of <num_primers> randomised primers with
    length <= <max_val>."""
    spacers = []
    for _ in range(num_primers):
        spacers.append((random.randrange(0, max_val),
                        random.randrange(0, max_val),
                        random.randrange(0, max_val),
                        random.randrange(0, max_val)))
    return spacers


class SeqFixtureManager:
    """Manages various datastructures required for testing."""
    incomplete_forward_primer: MBPrimerBuilder
    incomplete_reverse_primer: MBPrimerBuilder

    for_num_hetero: int
    rev_num_hetero: int

    for_max_spacer_length: int
    rev_max_spacer_length: int

    pot_forward_spacers: List[Tuple[int, int, int, int]]
    pot_reverse_spacers: List[Tuple[int, int, int, int]]

    forward_spacer: Tuple[int, int, int, int]
    reverse_spacer: Tuple[int, int, int, int]

    forward_spacer_seqs = []
    reverse_spacer_seqs = []

    num_to_generate: int

    def __init__(self) -> None:
        self.rsg4 = RandomSpacerGen(4, 4)
        self.rsg12 = RandomSpacerGen(12, 12)
        self.rsg_12_high_rigour = RandomSpacerGen(12, 12, 1)

        self.incomplete_forward_primer = MBPrimerBuilder()
        self.incomplete_reverse_primer = MBPrimerBuilder()

        self.for_num_hetero = 12
        self.rev_num_hetero = 12

        self.for_max_spacer_length = 12
        self.rev_max_spacer_length = 12

        self.pot_forward_spacers = []
        self.pot_reverse_spacers = []

        self.forward_spacer = (1, 1, 1, 1)
        self.reverse_spacer = (1, 1, 1, 1)

        self.forward_spacer_seqs = []
        self.reverse_spacer_seqs = []

        self.num_to_generate = LARGE_SAMPLE_SIZE

    def do_all(self):
        """Constructs all contained attributes."""
        self.set_primers_random()
        self.set_pot_spacers()
        self.set_single_spacers_random()
        self.set_spacers_seqs_random()

    def set_primers_random(self):
        """Sets the <incomplete_forward_primer> and <incomplete_reverse_primer>
        to random primers with region length generated as specified by defaults.
        """
        forward_len = random.randrange(BINDING_MIN, BINDING_MAX)
        reverse_len = random.randrange(BINDING_MIN, BINDING_MAX)
        self.incomplete_forward_primer = gen_incomplete_primer(forward_len)
        self.incomplete_reverse_primer = gen_incomplete_primer(reverse_len)

    def set_pot_spacers(self):
        """Sets <pot_forward_spacer> and <pot_reverse_spacer> to a set of
        heterogeneity spacer alignments."""
        sag = SpacerAlignmentGen(self.for_max_spacer_length,
                                 self.for_num_hetero)
        self.pot_forward_spacers = sag.get_all_spacer_combos(
            self.incomplete_forward_primer.get_binding_seq())

        sag.set_params(self.rev_max_spacer_length, self.rev_num_hetero)
        self.pot_reverse_spacers = sag.get_all_spacer_combos(
            self.incomplete_reverse_primer.get_binding_seq())

    def set_single_spacers_random(self):
        """Sets <self.forward_spacer> and <self.reverse_spacer> to a random
        spacer from the forward and reverse spacers."""
        self.forward_spacer = random.choice(self.pot_forward_spacers)
        self.reverse_spacer = random.choice(self.pot_reverse_spacers)

    def set_spacers_seqs_random(self):
        """Generates random forward and reverse spacers seqs using the
        attributes above."""
        gen_random_spacers(self.incomplete_forward_primer,
                           self.incomplete_reverse_primer,
                           self.forward_spacer_seqs,
                           self.reverse_spacer_seqs,
                           self.forward_spacer,
                           self.reverse_spacer,
                           self.num_to_generate)

# Begin fixtures for seq_alignment_analyser


def gen_seq_list(num_seqs: int, len: int, variance: int = 0,
                 degen: bool = False) -> List[Seq]:
    seqs = []
    for i in range(num_seqs):
        seq_len = random.randrange(int(len - variance/2), int(len + variance/2))
        if degen:
            seqs.append(gen_random_seq_degen(seq_len))
        else:
            seqs.append(gen_random_seq(seq_len))
    return seqs


T = TypeVar('T')


def all_unique(lst: List[T]) -> bool:
    """Where lst contains attributes with the __eq__ method defined."""
    for i_tup in enumerate(lst):
        i = i_tup[0]
        item_i = i_tup[1]
        for item_j in lst[i + 1:]:
            if item_j == item_i:
                return False
    return True
