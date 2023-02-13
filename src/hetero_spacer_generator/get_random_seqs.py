import random
from typing import Dict, List
from Bio.Seq import Seq

from src.hetero_spacer_generator.defaults import NUM_SPACERS
from src.hetero_spacer_generator.spacer_generator.spacer_filters import \
    SpacerAlignment, SpacerSet
from src.hetero_spacer_generator.primer_tools import IncompatibleSpacerError, MBPrimerBuilder, \
    PrimerSet
from src.hetero_spacer_generator.sequence_tools import remove_degen


def get_vacant_bases(spacers: SpacerAlignment) -> List[List[int]]:
    """Returns a List containing a List of indices of unfilled bases in
    <sequence_array>."""
    unfilled_bases = []
    for i in range(spacers[3] + 1):
        unfilled_bases.append([])

    for i in range(NUM_SPACERS):
        for j in range(spacers[i]):
            unfilled_bases[j].append(i)

    return unfilled_bases


def get_potential_bases(sequence_array: List[List[str]],
                        column: int) -> List[str]:
    """Returns  a list of all the bases that, if inserted into column, would
    maintain heterogeneity along the aligned nucleotide column."""
    possible_bases = ['A', 'T', 'C', 'G']
    for i in range(NUM_SPACERS):
        # The user entered a set of primers that do not guarantee
        # heterogeneity across the heterogeneity region.
        if sequence_array[i][column].isalnum() \
                and sequence_array[i][column] not in possible_bases and \
                sequence_array[i][column] != 'I':
            raise IncompatibleSpacerError()
        elif sequence_array[i][column].isalnum() and \
                sequence_array[i][column] != 'I':
            possible_bases.remove(sequence_array[i][column])
    return possible_bases


def gen_heterogeneity_spacers_rand(binding_seq: Seq,
                                   spacers: SpacerAlignment):
    """Returns a list of randomly generated valid heterogeneity spacers.
    When possible, avoids base repeats on same strand."""

    binding_seq = remove_degen(binding_seq)

    sequence_array = gen_sequence_array(binding_seq, spacers)
    unfilled_bases = get_vacant_bases(spacers)

    # Iterate through each column in reverse.
    for column in range(len(unfilled_bases) - 2, -1, -1):
        potential_bases = get_potential_bases(sequence_array, column)
        # Randomly fill one of the unspecified bases
        while unfilled_bases[column]:
            select_and_set(potential_bases, unfilled_bases, column,
                           sequence_array)

    # Extract newly generated heterogeneity spacers from the sequence array.
    final_spacers = []
    for i in range(NUM_SPACERS):
        spacer = Seq(str.join('', sequence_array[i][0:spacers[i]]))
        final_spacers.append(spacer)

    # Assumes NUM_SPACERS is 4
    return tuple(final_spacers)


def gen_hetero_set(incomplete_primer: MBPrimerBuilder,
                   spacers: SpacerAlignment, num_to_gen: int) \
        -> List[SpacerSet]:
    """Generates <num_to_gen> semi-random heterogeneity spaces.
    Ensures that primers match the lengths specified in <spacers> and ensure
    heterogeneity with the <incomplete_primer>'s binding_seq."""

    random_spacers = []
    binding_seq = incomplete_primer.get_binding_seq()
    for i in range(num_to_gen):
        random_spacers.append(
            gen_heterogeneity_spacers_rand(binding_seq, spacers))

    return random_spacers

def visualise_complete_primers(primer_sets: List[PrimerSet]) \
        -> Dict[int, PrimerSet]:
    """Prints all of <primers> to the presenter. Returns a dic mapping the
    value corresponding to each primer printed to that primer."""
    to_print = ''
    p_set_dict = {}
    for i in range(len(primer_sets)):
        p_set_dict[i + 1] = primer_sets[i]
        to_print += ''.join(["Primer Set #", str(i + 1), '\n'])
        to_print += str(primer_sets[i])
    return p_set_dict


def select_and_set(potential_bases: List[str],
                   unfilled_bases: List[List[int]],
                   column: int, sequence_array: List[List[str]]) -> None:
    """Selects a random base from <potential_bases> and places it in a
    random base position in <sequence_array>[base position][<column>]
    specified in <unfilled_bases>[<column>]. Removes the selected base from
    potential bases."""
    specific_potential_bases = potential_bases.copy()
    base_index = random.choice(unfilled_bases[column])

    # Avoid repeating the last base in the sequence if possible
    previous_base = sequence_array[base_index][column + 1]
    if len(specific_potential_bases) > 1 and \
            previous_base in specific_potential_bases:
        specific_potential_bases.remove(previous_base)

    selected_base = random.choice(specific_potential_bases)
    sequence_array[base_index][column] = selected_base
    unfilled_bases[column].remove(base_index)
    potential_bases.remove(selected_base)


def gen_sequence_array(binding_seq: Seq,
                       spacers: SpacerAlignment) \
        -> List[List[str]]:
    """Returns a representation of the <spacers> and <binding_seq>."""
    sequence_array = []  # [row][column]
    for i in range(len(spacers)):
        sequence_array.append([])
        # Fill hetero spacer region with empty slots
        for _ in range(spacers[i]):
            sequence_array[i].append('')  # Indicates an empty slot

        # Occupy remaining region (plus an extra base) with primer
        # Third spacer is always the largest.
        for j in range(0, spacers[3] - spacers[i] + 1):
            if j < len(binding_seq):
                sequence_array[i].append(binding_seq[j])
            else:
                sequence_array[i].append('-')

    return sequence_array
