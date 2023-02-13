from src.hetero_spacer_generator.primer_tools import MBPrimerBuilder
from typing import Union, Tuple
from src.hetero_spacer_generator.defaults import ABSOLUTE_MAX_NUM_HETERO, \
    ABSOLUTE_MAX_SPACER_LENGTH
from src.hetero_spacer_generator.spacer_generator.hetero_spacer_generator import HeteroGen

STR = 'str'
RANGE = 'ran'
DNA = 'dna'
DNA_DGN = 'dna_degen'
ACCEPT_DEGEN = False

def get_primer_and_spacers(hg: HeteroGen, direction: str,
                           is_one_step: bool = False) \
        -> Tuple[MBPrimerBuilder, Tuple[int, int, int, int]]:
    succ = False
    spacers = []
    incomplete_primer = MBPrimerBuilder()
    while not succ:
        print(
            "Please enter below the components of your " + direction + " primer.\n"
                                                                       "The adapter and binding regions are mandatory.\n")

        incomplete_primer = get_incomplete_primer(not is_one_step)

        print(
            "Please input the desired parameters for the construction of your "
            + direction + " primers.")

        spacer_length, num_hetero = get_params()
        hg.set_params(spacer_length, num_hetero)

        spacers = hg.get_all_spacer_combos(
            incomplete_primer.get_binding_seq())[0:50]

        if not spacers:
            print(
                "Failed to find a way in which to align the given sequences. Try"
                " setting the max spacer length to a greater value.")
        else:
            succ = True

    print("Spacers found. Press enter to select desired spacer combo.\n")

    num_to_spacer = hg.visualise_spacer_alignments(spacers,
                                                   incomplete_primer.
                                                   get_binding_seq())

    max_val = max(num_to_spacer.keys())

    spacer_ind = while_not_valid(
        "Enter the number of the primer alignment you'd "
        "like to use:",
        "Invalid input: Ensure the number you've entered "
        "represents one of the primers displayed.\n "
        "Enter the number of the primer alignment you'd "
        "like to use:", RANGE, start=1, end=max_val)

    spacer = num_to_spacer[int(spacer_ind)]
    return incomplete_primer, spacer


def get_spacer(hg: HeteroGen, incomplete_primer: MBPrimerBuilder,
               auto: bool = False) -> Tuple[int, int, int, int]:
    """Gets a spacer given some primer by showing a menu to the user."""
    spacers = hg.get_all_spacer_combos(
        incomplete_primer.get_binding_seq())[0:15]

    if auto:
        return spacers[0]

    if not spacers:
        raise ValueError(
            "Failed to find a way in which to align the given sequences. Try"
            " setting the max spacer length to a greater value.")

    print("Spacers found. Press enter to select desired spacer combo.\n")

    num_to_spacer = hg.visualise_spacer_alignments(spacers,
                                                   incomplete_primer.
                                                   get_binding_seq())

    max_val = max(num_to_spacer.keys())

    spacer_ind = while_not_valid(
        "Select the sizes of the heterogeneity spacers you'd like to use:",
        "Invalid input: Ensure the number you've entered "
        "represents one of the primers displayed.\n "
        "Enter the number of the primer alignment you'd "
        "like to use:", RANGE, start=1, end=max_val)

    return num_to_spacer[int(spacer_ind)]


def valid_input(input: Union[int, str], mode: str, start: Union[int, str] = 0,
                end: Union[int, str] = 0,
                allowed: Union[int, str] = '') -> bool:
    """Returns whether the input is an allowed input according to some
    parameters.
    Modes:
    STR:
        <input> is valid iff <input> in <allowed>
    RANGE:
        <input> is valid iff <start> =< <input> <= <end>.
        all of <input>, <start>, and <end> must be integers or string integers.
        """
    if mode == RANGE:
        start = int(start)
        end = int(end)
        neg_ind = input.find('-')
        if neg_ind >= 0:
            return input[1:].isnumeric() and start <= int(input) <= end
        return input.isnumeric() and start <= int(input) <= end
    if mode == STR:
        return str(input) in str(allowed)
    if mode == DNA:
        for char in input:
            if char.capitalize() not in 'ATCG':
                return False
        return True
    if mode == DNA_DGN:
        for char in input:
            if char.capitalize() not in 'ACTGUWSMKRYBDHVNI':
                return False
        return True


def while_not_valid(msg: str, e_msg: str, mode: Union[int, str],
                    start: Union[int, str] = '0',
                    end: Union[int, str] = '0',
                    allowed: Union[int, str] = '') -> str:
    """Continually prompts the user to enter input with the given <e_msg> while
    valid_input run with the given parameters is not true. Will try once with
    <msg> before printing <e_msg>."""
    inpt = input(msg)
    while not valid_input(inpt, mode, start, end, allowed):
        inpt = input(e_msg)
    return inpt


def get_incomplete_primer(get_index: bool = True) -> MBPrimerBuilder:
    """Prompts the user through the console for the components required to build
    a metabarcoding primer."""
    if ACCEPT_DEGEN:
        accepted = DNA_DGN
    else:
        accepted = DNA
    primer = MBPrimerBuilder()
    inpt = while_not_valid(
        "Enter the adapter region of your primer (5' - 3'): ",
        "Bad input, please ensure your entry is a valid sequence:", accepted)
    primer.set_adapter_seq(inpt.upper())

    if get_index:
        inpt = while_not_valid(
            "Enter the indexing region of your primer (5' - 3'): ",
            "Bad input, please ensure your entry is a valid sequence:",
            accepted)
        primer.set_index_seq(inpt.upper())

    inpt = while_not_valid(
        "Enter the binding region of your primer (5' - 3'): ",
        "Bad input, please ensure your entry is a valid sequence:", accepted)
    primer.set_binding_seq(inpt.upper())

    return primer


def get_params() -> Tuple[int, int]:
    """Prompts the user through the console for the parameters required to build
    a metabarcoding primer.
    Returns: Tuple[max spacer length, heterogeneity region size]"""
    spacer_length = while_not_valid(
        "Enter the maximum size of any heterogeneity spacer to be generated: ",
        "Bad input, please ensure your entry is a number greater than 0:",
        RANGE, start=0, end=ABSOLUTE_MAX_SPACER_LENGTH)

    num_hetero = while_not_valid(
        "Enter the size of the heterogeneity region you'd like in your primer: ",
        "Bad input, please ensure your entry is a number greater than 0:",
        RANGE, start=0, end=ABSOLUTE_MAX_NUM_HETERO)

    return int(spacer_length), int(num_hetero)

