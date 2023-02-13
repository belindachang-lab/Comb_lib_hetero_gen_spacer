
from typing import List, Dict, Callable, Union
from Bio.Seq import Seq
from src.hetero_spacer_generator.primer_tools import HeteroSeqTool, co_sort
from src.hetero_spacer_generator.sequence_tools import compare_bases_can_ignore, \
    remove_degen
from src.hetero_spacer_generator.spacer_generator.random_spacer_generator import RandomSpacerGen
from src.hetero_spacer_generator.primer_tools import MBPrimerBuilder, PrimerSet
from src.hetero_spacer_generator.spacer_generator.spacer_filters import \
    SpacerAlignment
from src.hetero_spacer_generator.defaults import GET_UNIQUE_SPACER_OFFSETS_DEFAULT, \
    NUM_SPACERS, MAX_SPACER_LENGTH, NUM_HETERO, \
    GET_SMALLEST_TOTAL_LEN_DEFAULT, GET_SMALLEST_OF_ANY_SPACER_DEFAULT, RIGOUR


# Begin methods for sorting spacer alignments
def get_smallest_total_len_list(spacers: List[SpacerAlignment]) -> List[int]:
    """Returns a list mapping the index of a spacer combo to the combined
    length of all spacers in that combo. """
    spacer_combined_length = []
    for spacer in spacers:
        spacer_combined_length.append(sum(spacer))
    return spacer_combined_length


def get_unique_spacer_offsets(spacers: List[SpacerAlignment]) -> List[int]:
    """Returns a list containing the number of different spacer lengths in
    each alignment."""
    num_unique = []
    for spacer in spacers:
        num_unique.append(4 - len(set(spacer)))
    return num_unique


def get_smallest_of_any_spacer_list(spacers: List[SpacerAlignment]) \
        -> List[int]:
    """Returns a list mapping the index of a spacer combo to the maximum length
    of a spacer in that combo."""
    spacer_max_len = []
    for spacer in spacers:
        spacer_max_len.append(max(spacer))
    return spacer_max_len


class SpacerAlignmentGen(HeteroSeqTool):
    """_max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across.
    _presenter:
            The class that handles the visualisation of data produced.
    _criteria:
            Functions score each spacer combo and return a list of those scores.
            Lower scores are always better.
    _criterion_to_weight:
            A dictionary mapping a criterion to it's relative weight (how
            important the scores returned by that function are when sorting
            spacer combos)."""

    _max_spacer_length: int
    _num_hetero: int
    _criteria: List[Callable[[List[SpacerAlignment]], List[int]]]
    _criterion_to_weight: \
        Dict[Callable[[List[SpacerAlignment]], List[int]], int]

    def _build_criteria(self) -> None:
        """Sets self._criteria to some list of default criteria and
        self._criterion_to_weight to the default weights for those criteria."""

        # Add new criteria here.
        self._criterion_to_weight = {
            get_smallest_of_any_spacer_list: GET_SMALLEST_OF_ANY_SPACER_DEFAULT,
            get_smallest_total_len_list: GET_SMALLEST_TOTAL_LEN_DEFAULT,
            get_unique_spacer_offsets: GET_UNIQUE_SPACER_OFFSETS_DEFAULT
        }

        self._criteria = list(self._criterion_to_weight.keys())

    def __init__(self, max_spacer_length: int, num_hetero: int) -> None:
        """Initialises the attributes to the values specified.
        Precondition:
            max_spacer_length >= num_hetero"""
        super().__init__(max_spacer_length, num_hetero)
        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero
        self._build_criteria()

    def get_all_spacer_combos(self, seq: Union[Seq, str]) \
            -> List[SpacerAlignment]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths. Spacer tuples are ordered from smallest spacer to
        greatest."""
        if not isinstance(seq, Seq):
            seq = Seq(seq)
        valid_spacer_combos = []
        spacers = []
        seq = remove_degen(seq)
        self._get_all_compatible_spacers(seq, spacers, NUM_SPACERS,
                                         spacer_combo_list=valid_spacer_combos)
        return valid_spacer_combos

    def _get_all_compatible_spacers(self, seq: Seq, spacers: List[int],
                                    target_length: int,
                                    spacer_combo_list: List[
                                        SpacerAlignment] = None) \
            -> None:
        """Will find all valid spacer combos for <seq> that include <spacers>.
        Will add all valid spacers with a length equal to <target_length> to
        <spacer_combo_list>."""

        if len(spacers) == target_length:
            if spacer_combo_list is not None:
                spacer_combo_list.append(tuple(spacers))
            return None

        for spacer in self._get_compatible_spacers(seq, spacers):
            # The set of spacers to be tested for heterogeneity across
            # self.num_hetero
            trial_spacers = spacers.copy()
            trial_spacers.append(spacer)
            self._get_all_compatible_spacers(seq, trial_spacers, target_length,
                                             spacer_combo_list)
        return None

    def _get_compatible_spacers(self, seq: Seq,
                                seqs_spacers: List[int]) -> tuple[int, ...]:
        """Returns a tuple containing all spacers < <self.num_hetero> such
        that for all j + spacer < <self.num_hetero>, seq[j + spacer] != any
        seqs[i][j]. Returns an empty tuple if no such spacers exist."""

        valid_spacers = []
        # Assume that there exist no valid spacers less than the greatest
        # already used
        if len(seqs_spacers) > 0:
            start_point = seqs_spacers[-1]
        else:
            start_point = 0

        # Iterate though all possible spacers - add to valid spacers if it is
        # compatible
        for i in range(start_point, self._max_spacer_length + 1):
            if self._is_compatible_spacer(seq, seqs_spacers, i):
                valid_spacers.append(i)

        return tuple(valid_spacers)

    def _is_compatible_spacer(self, seq: Seq, spacers: List[int],
                              spacer: int, ) -> bool:
        """Returns whether <seq> shifted by <spacer> has any matching bases
        in the same position as any sequence in <seqs> - each shifted by their
        respective spacers.

        Preconditions: Every item in spacers is less than spacer.
        """
        # Checking for uniqueness at all bases
        for i in range(0, self._num_hetero - spacer):

            # Checking each of the seqs for matching bases.
            for j in range(0, len(spacers)):
                try:
                    original_spacer_base = seq[i + spacer - spacers[j]]
                    new_spacer_seq_base = seq[i]
                    if compare_bases_can_ignore(original_spacer_base,
                                                new_spacer_seq_base):
                        return False
                    # Sequence is shorter than heterogeneity region, spacer
                    # fails to push the sequence out of the heterogeneity
                    # region.
                except IndexError:
                    return False

        return True

    def sort_spacer_combos(self,
                           spacer_combos: List[SpacerAlignment]
                           ) -> None:
        """Returns a list of spacer combos sorted according to <self.criteria>.
        """
        # Instantiate a list of scores with the lowest possible scores. This
        # will store the scores for each spacer in <spacer_combos>
        scores = []
        for i in range(len(spacer_combos)):
            scores.append(1)

        for i in range(len(self._criteria)):
            criterion = self._criteria[i]
            weight = self._criterion_to_weight[criterion]
            criterion_scores = criterion(spacer_combos)
            for j in range(len(spacer_combos)):
                scores[j] = scores[j] * criterion_scores[j] * weight

        co_sort(scores, spacer_combos)

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero = num_hetero
        self._max_spacer_length = max_spacer_length


class HeteroGen:
    """Manages the the creation of heterogeneity spacers. Facades
    SpacerAlignmentGen and RandomPrimerGen.

    === Private Attributes ===
    _max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across.
    _presenter:
            The class that handles the visualisation of data produced.
    _alignment_gen:
            Produces alignments of the binding sequences so that heterogeneity
            can be ensured between them.
    _primer_gen:
            Generates primer sequences.
            """
    _max_spacer_length: int
    _num_hetero: int
    _alignment_gen: SpacerAlignmentGen
    _primer_gen: RandomSpacerGen

    def __init__(self, max_spacer_length: int = MAX_SPACER_LENGTH,
                 num_hetero: int = NUM_HETERO, rigour: int = RIGOUR) -> None:
        """Initialises the attributes to the values specified."""
        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero
        self._alignment_gen = SpacerAlignmentGen(max_spacer_length, num_hetero,)
        self._primer_gen = RandomSpacerGen(max_spacer_length, num_hetero,
                                           rigour=rigour)

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero = num_hetero
        self._max_spacer_length = max_spacer_length
        self._primer_gen.set_params(max_spacer_length, num_hetero)
        self._alignment_gen.set_params(max_spacer_length, num_hetero)

    def get_primer_gen(self) -> RandomSpacerGen:
        """Returns this HeteroSpacerGen's primer gen"""
        return self._primer_gen

    def set_rigour(self, rigour: int) -> None:
        """Sets the rigour to the specified value."""
        self._primer_gen.set_rigour(rigour)

    def set_pairwise(self, degen: bool = None) -> None:
        """Sets this HeteroGen to optimise the generated primers for pairwise
        PCR."""
        self._primer_gen.set_pairwise(degen=degen)

    def get_all_spacer_combos(self, seq: Seq
                              ) -> List[SpacerAlignment]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths. Will sort the alignments in order to provide shorter
        spacers."""
        aligns = self._alignment_gen.get_all_spacer_combos(seq)
        self._alignment_gen.sort_spacer_combos(aligns)
        return aligns

    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacer: SpacerAlignment,
                        reverse_spacer: SpacerAlignment,
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences."""
        return self._primer_gen.get_hetero_seqs(incomplete_forward_primer,
                                                incomplete_reverse_primer,
                                                forward_spacer, reverse_spacer,
                                                num_to_return)

    def visualise_spacer_alignments(self,
                                    spacers: List[SpacerAlignment],
                                    seq: Seq) -> Dict[int, SpacerAlignment]:
        """Displays a visual representation of all of the <spacers> for the
        given <seq>. Assumes a set of 4 spacers per tuple."""
        spacer_dict = {}
        to_print = ''
        for i in range(0, len(spacers)):
            # Add title
            spacer_tup = spacers[i]
            spacer_dict[i + 1] = spacers[i]
            to_print += 'Spacer #' + str(i + 1) + ' ' + str(spacer_tup) + '\n'

            # Add indexing numbers
            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            # Add sequences and spacers
            for j in range(0, NUM_SPACERS):
                spacer = spacer_tup[j]

                to_print += '+ ' * min(spacer, self._num_hetero)
                to_print += '+' * max(spacer - self._num_hetero, 0)

                if spacer < self._num_hetero:
                    to_print += ' '.join(
                        str(seq[0:self._num_hetero - spacer])) + ' '

                to_print += ''.join([''.join(
                    str(seq[self._num_hetero - spacer:])),
                    '\n'])

            to_print += '\n'

        print(to_print)
        return spacer_dict

    def visualise_seq_arr(self, seq_arrs: List[List[List[str]]]) -> None:
        """Displays a visual representation of a given <seq_arrs>.

        Precondition: For all valid i, j, and k,
        len(<seq_arr>[i][j][k]) > self._num_hetero

        Note:
            Replaces all '' in seq_arr with 'N' """
        to_print = ''
        for i in range(len(seq_arrs)):
            for j in range(len(seq_arrs[i])):
                for k in range(len(seq_arrs[i][j])):
                    if seq_arrs[i][j][k] == '':
                        seq_arrs[i][j][k] = 'N'

        for i in range(len(seq_arrs)):
            to_print += ''.join(["Sequence array #", str(i + 1), '\n'])

            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            # Add spaces between the bases in the heterogeneity region
            for j in range(NUM_SPACERS):
                to_print += ' '.join(seq_arrs[i][j][0:self._num_hetero]) + ' '
                to_print += ''.join(seq_arrs[i][j][self._num_hetero:]) + '\n'

            to_print += '\n'

        print(to_print)
