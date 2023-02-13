import functools
from abc import ABC, abstractmethod
import random
from statistics import stdev
from typing import Any, Callable, Collection, Iterator, List, Tuple, \
    Dict, Iterable, Union, TypeVar

from infinity import Infinity
from multiprocessing import Process
from src.hetero_spacer_generator.sequence_tools import get_max_complementarity, \
    get_max_complementarity_consec
from Bio.Seq import Seq

FORWARD = "forward"
REVERSE = "reverse"
TAB = '    '

DEFAULT_PS_TYPE = 'pai'
PAIRWISE = 'pai'
SIMULTANEOUS = 'sim'

# Lower values make variance a more important variable in score calculation
# decisions
VARIANCE_IMPORTANCE = 3


class Primer(Seq):
    """Basic methods common among all immutable primers"""

    def __init__(self, seq: str):
        super().__init__(seq)


class HeteroSeqTool(ABC):
    """A class that handles heterogeneity spacers.

    --- Private Attributes ---
    _max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across."""

    def __init__(self, max_spacer_length: int, num_hetero: int) -> None:
        """Initialises the attributes to the values specified.
        Precondition:
            max_spacer_length >= num_hetero"""
        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero


class MBPrimer(Primer):
    """Stores the attributes of a primer used in two step PCR for metabarcoding.

    === Private Attributes ===
    _adapter_seq:
            The adapter region of the primer.
    _index_seq:
            The indexing region.
    _heterogen_seq:
            The heterogeneity sequence of the primer.
    _binding_seq:
            The specific binding region.
    _direction:
            Primers are either "forward" or "reverse".
    _3p_len:
            The length of the region downstream of the heterogeneity spacer.
    _5p_len:
            The length of the region upstream of the heterogeneity spacer.

    On Directionality:
            The sequence of the primer will be given in the 5'-3' direction.
    """
    _adapter_seq: Seq
    _index_seq: Seq
    _heterogen_seq: Seq
    _binding_seq: Seq
    _direction: str
    _3p_len: int
    _5p_len: int

    def __init__(self, adapter_seq: str, index_seq: str, heterogen_seq: str,
                 binding_seq: str) -> None:
        """Constructs a primer sequence with components in the following order:
        <adapter_seq> <index_seq> <heterogen_seq> <primer_region> with each of
        their sequences read left to right."""
        super().__init__(''.join([adapter_seq, index_seq, heterogen_seq,
                                  binding_seq]))

        self._adapter_seq = Seq(adapter_seq)
        self._index_seq = Seq(index_seq)
        self._heterogen_seq = Seq(heterogen_seq)
        self._binding_seq = Seq(binding_seq)
        self._5p_len = len(adapter_seq) + len(index_seq)
        self._3p_len = len(binding_seq)

    def get_adapter_seq(self) -> Seq:
        """Returns adapter_seq."""
        return self._adapter_seq

    def get_index_seq(self) -> Seq:
        """Returns index_seq."""
        return self._index_seq

    def get_heterogen_seq(self) -> Seq:
        """Returns heterogen_seq."""
        return self._heterogen_seq

    def get_binding_seq(self) -> Seq:
        """Returns binding_seq."""
        return self._binding_seq

    def get_heterogen_seq_ind(self) -> int:
        """Returns the index of the first base of the heterogeneity sequence."""
        return len(self._adapter_seq)

    def get_5p(self) -> Seq:
        """Returns the sequence to the 5' end of the heterogeneity spacer."""
        return Seq(''.join([self._adapter_seq.__str__(),
                            self._index_seq.__str__()]))

    def get_3p(self) -> Seq:
        """Returns the sequence to the 3' end of the heterogeneity spacer."""
        return self._binding_seq

    def get_5p_len(self) -> int:
        """Returns the sequence to the 5' end of the heterogeneity spacer."""
        return self._5p_len

    def get_3p_len(self) -> int:
        """Returns the sequence to the 3' end of the heterogeneity spacer."""
        return self._3p_len


class MBPrimerBuilder(Iterable):
    """Stores the attributes of a primer used in two step PCR for metabarcoding.

    === Private Attributes ===
    _adapter_seq:
            The adapter region of the primer.
    _index_seq:
            The indexing region.
    _heterogen_seq:
            The heterogeneity sequence of the primer.
    _binding_seq:
            The locus specific binding region.
    _direction:
            Primers are either "forward" or "reverse".

    On Directionality:
            The sequence of the primer will be given in the 5'-3' direction.
    """

    _adapter_seq: Seq
    _index_seq: Seq
    _heterogen_seq: Seq
    _binding_seq: Seq
    _inherent_consec: int
    _inherent_total: int

    def __init__(self, adapter_seq: Seq = Seq(''), index_seq: Seq = Seq(''),
                 heterogen_seq: Seq = Seq(''), binding_seq: Seq = Seq('')) \
            -> None:
        """Initialises this MBPrimerBuilder with empty seqs unless specified."""
        self._adapter_seq = adapter_seq
        self._index_seq = index_seq
        self._heterogen_seq = heterogen_seq
        self._binding_seq = binding_seq

    def __iter__(self) -> Iterator[Seq]:
        return [self._adapter_seq, self._index_seq, self._heterogen_seq,
                self._binding_seq].__iter__()

    def get_seq(self) -> Seq:
        """Returns a seq representation of this primer."""
        return Seq(''.join([self._adapter_seq[:], self._index_seq[:],
                            self._heterogen_seq[:], self._binding_seq[:]]))

    def set_adapter_seq(self, seq: Union[Seq, str]) -> None:
        """Sets adapter_seq to <seq>."""
        self._adapter_seq = Seq(str(seq))

    def set_index_seq(self, seq: Union[Seq, str]) -> None:
        """Sets index_seq to <seq>."""
        self._index_seq = Seq(str(seq))

    def set_heterogen_seq(self, seq: Union[Seq, str]) -> None:
        """Sets heterogen_seq to <seq>."""
        self._heterogen_seq = Seq(str(seq))

    def get_heterogen_seq_ind(self) -> int:
        """Returns the index of the first base of the heterogeneity sequence."""
        return len(self._adapter_seq)

    def set_binding_seq(self, seq: Union[Seq, str]) -> None:
        """Sets binding_seq to <seq>."""
        self._binding_seq = Seq(str(seq))

    def get_adapter_seq(self) -> Seq:
        """Returns adapter_seq."""
        return self._adapter_seq

    def get_index_seq(self) -> Seq:
        """Returns index_seq."""
        return self._index_seq

    def get_heterogen_seq(self) -> Seq:
        """Returns heterogen_seq."""
        return self._heterogen_seq

    def get_binding_seq(self) -> Seq:
        """Returns binding_seq."""
        return self._binding_seq

    def get_5p(self) -> Seq:
        """Returns the sequence to the 5' end of the heterogeneity spacer."""
        return Seq(''.join([self._adapter_seq.__str__(),
                            self._index_seq.__str__()]))

    def get_3p(self) -> Seq:
        """Returns the sequence to the 3' end of the heterogeneity spacer."""
        return self._binding_seq

    def get_mbprimer(self) -> MBPrimer:
        """Returns a completed version of this primer."""
        return MBPrimer(self._adapter_seq.__str__(),
                        self._index_seq.__str__(),
                        self._heterogen_seq.__str__(),
                        self._binding_seq.__str__())


SimpleCriterion = Callable[[MBPrimer, MBPrimer], int]


def spacers_to_primers(incomplete_primer: MBPrimerBuilder,
                       spacers: Tuple[Seq, Seq, Seq, Seq]) -> List[MBPrimer]:
    """Completes the <incomplete_primer> with each of <spacers> and returns  a
    list of the results."""
    primers = []
    for spacer in spacers:
        incomplete_primer.set_heterogen_seq(spacer)
        primers.append(incomplete_primer.get_mbprimer())
    return primers


def remove_highest_scores(lst: List[Any], scores_dict: Dict[int, List[int]],
                          num_to_keep: int, lowest: bool = False) -> None:
    """Given <scores_dict> which maps scores to a list of indices, removes items
    at all indices in <lst> except for the <num_to_keep> highest scoring items.
    Iff <lowest> is true, then this behaviour is reversed, and the lowest
    scoring items will be kept."""
    best_to_worst = []
    sorted_keys = list(scores_dict.keys())
    sorted_keys.sort()
    # Insert indices such that best_to_worse contains indices in the
    # lowest scoring -> highest scoring direction
    for key in sorted_keys:
        for index in scores_dict[key]:
            best_to_worst.append(index)

    # Score preference is reversed if lower scores are preferred.
    if lowest:
        best_to_worst.reverse()

    # Remove unwanted items at higher indices first.
    ind_to_remove = best_to_worst[num_to_keep: len(best_to_worst)]
    ind_to_remove.sort(reverse=True)

    for index in ind_to_remove:
        lst.pop(index)


class HalfSet:
    """A set of 4 forward or reverse primers.

    === Public Attributes ===
    primers:
            A set of valid Metabarcoding primers"""

    primers: Tuple[MBPrimer, MBPrimer, MBPrimer, MBPrimer]
    _average_score: float
    _num_scores: int
    _is_active: bool

    def __init__(self, incomplete_primer: MBPrimerBuilder,
                 spacers: Tuple[Seq, Seq, Seq, Seq]) -> None:
        """Constructs a set of primers complete using <incomplete_primer> each
        with one of <spacers>"""
        self._average_score = 0
        self._num_scores = 0
        self._is_active = True
        primers = []
        for spacer in spacers:
            incomplete_primer.set_heterogen_seq(spacer)
            primers.append(incomplete_primer.get_mbprimer())
        incomplete_primer.set_heterogen_seq(Seq(''))
        self.primers = (primers[0], primers[1], primers[2], primers[3])

    def __iter__(self):
        """Returns an iterator for the set of primers contained in this set."""
        return iter(self.primers)

    def deactivate(self) -> None:
        """Deactivates this HalfSet."""
        self._is_active = False

    def is_active(self) -> bool:
        """Returns whether this halfset is still active."""
        return self._is_active

    def get_avg(self) -> float:
        """Returns the average score of this set."""
        return self._average_score

    def set_avg(self, new_avg: float) -> None:
        """Sets the average to the given value"""
        self._average_score = new_avg
        return

    def reset_avg(self) -> None:
        """Resets the values associated with calculating this HalfSets avg
        performance."""
        self._average_score = 0
        self._num_scores = 0

    def update_avg_score(self, new_score: Union[int, List[int]]) -> None:
        """Adds <new_score> to this HalfSets running average"""
        if type(new_score) == int:
            self._num_scores += 1
            self._average_score = self._average_score * \
                                  (self._num_scores - 1) / self._num_scores + \
                                  new_score / self._num_scores
        elif type(new_score) == list:
            old_score = self._num_scores
            self._num_scores += len(new_score)
            self._average_score = self._average_score * \
                                  old_score / self._num_scores + \
                                  new_score / self._num_scores


class PrimerSet(Iterable):
    """A set of MBPrimers.
    === Private Attributes ===
    _forward_primers:
            The forward primers contained in this set.
    _reverse_primers:
            The reverse primers contained by this set."""

    _forward_primers: List[MBPrimer]
    _reverse_primers: List[MBPrimer]

    def __init__(self, forward_primers: Iterable[MBPrimer] = (),
                 reverse_primers: Iterable[MBPrimer] = ()) -> None:
        """Initialises the Primer set with the given list of primers."""
        self._forward_primers = []
        self._reverse_primers = []
        for primer in forward_primers:
            self._forward_primers.append(primer)
        for primer in reverse_primers:
            self._reverse_primers.append(primer)

    def __iter__(self) -> List[Tuple[MBPrimer, MBPrimer]]:
        """Returns an iterable list of this PrimerSets possible forward and
        reverse primer combinations."""
        lst = []
        for forward in self._forward_primers:
            for reverse in self._reverse_primers:
                lst.append((forward, reverse))
        return lst

    def get_forward_primers(self) -> List[MBPrimer]:
        """Returns this sets forward primers."""
        return self._forward_primers

    def get_reverse_primers(self) -> List[MBPrimer]:
        """Returns this sets forward primers."""
        return self._forward_primers

    def get_score(self) -> float:
        return -1

    # TODO THIS IS BAD FIX THIS
    def apply_criteria(self, criteria: List[SimpleCriterion],
                       weights: List[int]) -> int:
        pass

    def get_plain_seqs(self) -> str:
        """Returns the sequences of this primer set as plain strings seperated
        by the newline character."""
        str_rep = ''
        for primer in self._forward_primers:
            str_rep += ''.join([str(primer.get_adapter_seq()),
                                str(primer.get_index_seq()),
                                str(primer.get_heterogen_seq()),
                                str(primer.get_binding_seq())])
            str_rep += '\n'
        for primer in self._reverse_primers:
            str_rep += ''.join([str(primer.get_adapter_seq()),
                                str(primer.get_index_seq()),
                                str(primer.get_heterogen_seq()),
                                str(primer.get_binding_seq())])
            str_rep += '\n'
        return str_rep

    def get_fasta_seqs(self, set_number: int = 0) -> str:
        """Returns the sequences of this primer set as plain strings seperated
        by the newline character. Will name each primer according to the
        <set_number>, where each set is assumed to contain 4 primers. Does start
        at 0."""
        primer_number = set_number * 4 + 1
        str_rep = ''
        for primer_inds in range(len(self._forward_primers)):
            for_primer = self._forward_primers[primer_inds]
            rev_primer = self._reverse_primers[primer_inds]
            str_rep += ''.join(['>Forward #', str(primer_number), '\n'])
            str_rep += ''.join([str(for_primer.get_adapter_seq()),
                                str(for_primer.get_index_seq()),
                                str(for_primer.get_heterogen_seq()).lower(),
                                str(for_primer.get_binding_seq())])
            str_rep += '\n'
            str_rep += ''.join(['>Reverse #', str(primer_number), '\n'])
            str_rep += ''.join([str(rev_primer.get_adapter_seq()),
                                str(rev_primer.get_index_seq()),
                                str(rev_primer.get_heterogen_seq()).lower(),
                                str(rev_primer.get_binding_seq())])
            str_rep += '\n'
            primer_number += 1
        return str_rep

    def __str__(self) -> str:
        """Returns a string representation of this primer set."""
        str_rep = ''
        str_rep += 'Forward Primers:\n'
        for primer in self._forward_primers:
            str_rep += TAB + TAB.join([str(primer.get_adapter_seq()),
                                       str(primer.get_index_seq()),
                                       str(primer.get_heterogen_seq()),
                                       str(primer.get_binding_seq())])
            str_rep += '\n'
        str_rep += 'Reverse Primers:\n'
        for primer in self._reverse_primers:
            str_rep += TAB + TAB.join([str(primer.get_adapter_seq()),
                                       str(primer.get_index_seq()),
                                       str(primer.get_heterogen_seq()),
                                       str(primer.get_binding_seq())])
            str_rep += '\n'
        return str_rep


SpacerPairing = Tuple[int, int, int, int]


def calculate_score(scores: Collection[int]) -> int:
    """Calculates the average score of <scores>, increasing the score if
    <scores> has high variance."""
    std = stdev(scores)
    avg = sum(scores) / len(scores)
    return int(avg + std / VARIANCE_IMPORTANCE)


# Where <SpacersSets> is the set of all spacers seqs, <MBPrimerBuilder> is the
# incomplete spacer. Returns a list of all scores.
PairWiseCriterionSingle = Callable[[MBPrimer], int]


class PairwisePrimerSet(PrimerSet):
    """A set of MBPrimers intended to be used in pairs.
    --- Private Attributes ---
    _optimal_pairing:
            A pair of some of this sets forward and reverse primers that have
            minimal hetero-binding. Where the pairs are:
            forward_primers[i] pairs with reverse_primers[optimal_pairing[i]].
    _pairing_scores:
            The scores produced by the given pairings of primers.
    _min_pairing_score:
            The best score produced by any pairing.
    _has_been_scored:
            True iff this set has been evaluated by criteria.
    """
    _optimal_pairing: SpacerPairing
    _pairing_scores: Dict[SpacerPairing, int]
    _min_pairing_score: int
    _been_scored: bool

    def __init__(self, forward_primers: Iterable[MBPrimer] = (),
                 reverse_primers: Iterable[MBPrimer] = ()) -> None:
        """Initialises the Primer set with the given list of primers."""
        super().__init__(forward_primers, reverse_primers)
        self._optimal_pairing = (0, 1, 2, 3)
        self._pairing_scores = {}
        self._min_pairing_score = 0
        self._been_scored = False

        for pairing in POSSIBLE_PAIRINGS:
            self._pairing_scores[pairing] = 0

    def __str__(self) -> str:
        """Returns a string representation of this primer set."""
        str_rep = ''
        str_rep += super().__str__()
        str_rep += "Recommended Pairing:" + \
                   ', '.join(["F1 - R{R1:d}",
                              "F2 - R{R2:d}",
                              "F3 - R{R3:d}",
                              "F4 - R{R4:d}"]).format(
                       R1=self._optimal_pairing[0] + 1,
                       R2=self._optimal_pairing[1] + 1,
                       R3=self._optimal_pairing[2] + 1,
                       R4=self._optimal_pairing[3] + 1, )
        return str_rep

    def __iter__(self) -> List[Tuple[MBPrimer, MBPrimer]]:
        """Returns an iterable list of this PrimerSets forward and reverse
        primers."""
        lst = []
        for f, r in enumerate(self._optimal_pairing):
            lst.append((self._forward_primers[f], self._reverse_primers[r]))
        return lst

    def randomise_optimal_pairing(self) -> None:
        """Randomises this sets pairwise scores."""
        lst = [0, 1, 2, 3]
        random.shuffle(lst)
        pairing = tuple(lst)
        self._optimal_pairing = pairing
        self._min_pairing_score = self._pairing_scores[pairing]

    def get_score(self) -> float:
        return self._min_pairing_score

    def get_min_pairing_score(self) -> int:
        """Gets the minimum pairing score of this set."""
        return self._min_pairing_score

    def has_been_scored(self) -> bool:
        """Returns whether this set has been scored."""
        return self._been_scored

    def apply_criteria(self, criteria: List[SimpleCriterion],
                       weights: List[int]) -> int:
        """Calculates the scores for all parings using all of <criteria> and
        weights, returning the minimum of those scores.
        Assumes:
            Each criteria[i] corresponds to weights[i]
            i.e. len(criteria) == len(weights)."""
        for i in range(len(criteria)):
            self.apply_criterion(criteria[i], weights[i])
        self._been_scored = True
        return self.update_min_score()

    def update_min_score(self) -> int:
        """Updates <self.min_score> to the minimum score in pairing_scores."""
        self._min_pairing_score = min(self._pairing_scores.values())
        self._optimal_pairing = \
            {v: k for k, v in self._pairing_scores.items()}[
                self._min_pairing_score]
        return self._min_pairing_score

    def apply_criterion(self, criterion: SimpleCriterion, weight: int) -> None:
        """Updates <self.pairing_scores> using <criterion>."""
        for pairing in POSSIBLE_PAIRINGS:
            scores = []
            for i in range(4):
                scores.append(criterion(self._forward_primers[i],
                                        self._reverse_primers[pairing[i]]))
            self._pairing_scores[pairing] += calculate_score(scores) \
                                             * weight


def eval_total_complementarity(incomplete_primer: MBPrimerBuilder,
                               spacers: Tuple[Seq, Seq, Seq, Seq]) -> int:
    """Returns the max complementarity between any of <spacers> and any
    <incomplete_primer> completed with a spacer in <spacers>."""
    cmplmnt_lst = []
    primers = spacers_to_primers(incomplete_primer, spacers)
    for spacer in spacers:
        # Calculate max_comp for any forward primer and forward hetero.
        rev_hetero = spacer.__str__()[::-1]
        cmplmnt_lst.append(get_max_complementarity(Seq(rev_hetero),
                                                   primers))
        # No need to reverse hetero seq when calculating binding to reverse
        # primer.
    return max(cmplmnt_lst)


def eval_consecutive_complementarity(incomplete_primer: MBPrimerBuilder,
                                     spacers: Tuple[Seq, Seq, Seq, Seq]) -> int:
    """Returns the max complementarity between any of <spacers> and any
    <incomplete_primer> completed with a spacer in <spacers>."""
    cmplmnt_lst = []
    primers = spacers_to_primers(incomplete_primer, spacers)
    for spacer in spacers:
        # Calculate max_comp for any forward primer and forward hetero.
        rev_hetero = spacer.__str__()[::-1]
        cmplmnt_lst.append(get_max_complementarity_consec(Seq(rev_hetero),
                                                          primers))
        # No need to reverse hetero seq when calculating binding to reverse
        # primer.
    return max(cmplmnt_lst)


def evaluate_heterogen_binding_cross(forward_primers: HalfSet,
                                     reverse_primers: HalfSet) -> int:
    """Returns the max site complementarity between the heterogeneity
    sequences in <forward_primers> and the <reverse_primers> and visa versa. """
    max_comp = 0

    # Calculate the maximum complementarity between any of the reverse and
    # forward primers
    for primer in forward_primers:
        comp = get_max_complementarity(primer.get_heterogen_seq(),
                                       reverse_primers)
        if comp > max_comp:
            max_comp = comp

    # Same process for reverse primers
    for primer in reverse_primers:
        comp = get_max_complementarity(primer.get_heterogen_seq(),
                                       forward_primers)
        if comp > max_comp:
            max_comp = comp

    return max_comp


def co_sort(to_sort: List[int], to_follow: List[Any], reverse: bool = False) \
        -> None:
    """Sorts to_sort least to greatest (greatest to least if reverse is true).
    Whenever an item in to_sort is moved index i to index j, the item in
    to_follow at i is moved to j.

    Precondition:
            len(to_follow) >= len(to_sort)"""
    inv = [-1, 1][reverse]
    for i in range(len(to_sort), 1, -1):
        for j in range(i - 1):
            if to_sort[j] * inv < to_sort[j + 1] * inv:
                to_sort[j], to_sort[j + 1] = to_sort[j + 1], to_sort[j]
                to_follow[j], to_follow[j + 1] = to_follow[j + 1], to_follow[j]


def co_insert(values: List[int], to_follow: List[Any],
              val_to_insert: int, follow_to_insert: Any, reverse: bool = False) \
        -> None:
    """Inserts <val_to_insert> into a the sorted list (L-G) <values> such that
    <values> is still sorted. Will insert <follow_to_insert> into the same index
    in <to_follow> as <val_to_insert> was inserted into <values>."""
    # Used to reverse logic in value comparisons
    dir = 1 - 2 * reverse
    for i in range(len(values)):
        if val_to_insert * dir < values[i] * dir:
            values.insert(i, val_to_insert)
            values.pop()
            to_follow.insert(i, follow_to_insert)
            to_follow.pop()
            return


def get_arrangements_recursive(possible_states: List[int],
                               arrangements: List[Tuple[int]],
                               arrangement: List[int],
                               depth: int, target_depth: int) -> None:
    """Adds one of <possible_states> to <arrangement>, will recurse until
    <target_depth> is reached, adding an <arrangement> of length <target_depth>
    to <arrangements>."""
    if depth == target_depth:
        arrangements.append(tuple(arrangement))
        return

    # Recurse for each possible state
    for i in range(len(possible_states)):
        next_arrangement = arrangement + [possible_states[i]]
        next_possible_states = possible_states[:]
        next_possible_states.pop(i)
        next_depth = depth + 1
        get_arrangements_recursive(next_possible_states, arrangements,
                                   next_arrangement, next_depth, target_depth)
    return


def get_all_arrangements(num_states: int, num_positions: int):
    """Gets all possibles arrangements of a list of length <num_positions> with
    each value in the list being a unique number ST 0 =< n < num states
    Returns a List containing tuples of possible combinations where
    len(tuple) == num_positions"""

    possible_states = []
    arrangements = []
    for i in range(num_states):
        possible_states.append(i)
    get_arrangements_recursive(possible_states, arrangements, [], 0,
                               num_positions)
    return arrangements


POSSIBLE_PAIRINGS = get_all_arrangements(4, 4)


def add(dic: Dict[Any, List[Any]], key: Any, item: Any) -> None:
    """Adds <item> to the <dic> with <key> if <key> already maps to a list, else
     creates a list at <key> and adds <item> to it."""
    if key in dic.keys():
        dic[key].append(item)
    else:
        dic[key] = [item]


def get_n_lowest_matrix(scores: List[List[int]],
                        n: int, highest: bool = False) \
        -> Tuple[List[Tuple[int, int]], List[int]]:
    """Returns a tuple containing: a list of the indices of the <n> lowest
    scores in scores, and a list of those scores. Iff <highest>,
    will return the <n> highest scoring items.

    Assumes:
     - That every score has a correlated item.
     - An item's index in <items> also points to its score in <scores>"""

    # Best scores sorted from lowest to highest complementarity
    min_scores = []
    # The combination of for and rev primers that produced the scores in
    # min_scores.
    min_items = []
    # Set the scores to a list of the worst possible scores
    worst_score = MaxInt() if not highest else - MaxInt()

    rev = -2 * highest + 1

    for i in range(n):
        min_scores.append(worst_score)
        min_items.append((-1, -1))

    # If a score is better than the worst score in scores, replace that
    # score and combo with the better score and combo. Sort best to worst.
    for i in range(len(scores)):
        for j in range(len(scores[i])):
            # Iff highest reverse logic.
            if scores[i][j] * rev < min_scores[-1] * rev:
                co_insert(min_scores, min_items,
                          scores[i][j], (i, j), reverse=highest)

    return min_items, min_scores


def get_n_lowest(scores: List[Union[int, float]], n: int, highest: bool = False) \
        -> Tuple[List[int], List[int]]:
    """Returns a tuple containing [0]: a list of the indices of the <n> lowest
    scores in scores, [1]: a list of those scores. Iff <highest>,
    will return the <n> highest scoring items.

    Assumes:
     - That every score has a correlated item.
     - An item's index in <items> also points to its score in <scores>"""

    if n < 1:
        raise ValueError("n must be greater than 0")
    # Best scores sorted from lowest to highest complementarity
    min_scores = []
    # The combination of for and rev primers that produced the scores in
    # min_scores.
    min_items = []
    # Set the scores to a list of the worst possible scores
    worst_score = MaxInt() if not highest else - MaxInt()

    rev = -2 * highest + 1

    for i in range(n):
        min_scores.append(worst_score)
        min_items.append(-1)

    # If a score is better than the worst score in scores, replace that
    # score and combo with the better score and combo. Sort best to worst.
    for i in range(len(scores)):
        # Iff highest reverse logic.
        if scores[i] * rev < min_scores[-1] * rev:
            co_insert(min_scores, min_items,
                      scores[i], i, reverse=highest)

    return min_items, min_scores

def procs_join(procs: List[Process]) -> None:
    """Pauses the program until all of <procs> have joined."""
    for proc in procs:
        proc.join()

@functools.total_ordering
class Comparable(ABC):
    @abstractmethod
    def __lt__(self, other) -> bool:
        pass

    @abstractmethod
    def __eq__(self, other) -> bool:
        pass


T = TypeVar('T', Comparable, int, float)


def get_n_highest_sbs(vals: List[List[T]], n: int,
                      least: bool = False) -> List[T]:
    """Returns the <n> greatest Comparable[T]s in vals.

    Precondition:
        vals must be sorted G-L, (L-G iff <least>)
        combined length of all <vals> must be greater than or equal to <n>."""
    t_len = 0
    for val in vals:
        t_len += len(val)
    if t_len < n:
        raise ValueError("Combined length of all <vals> must be greater than or"
                         " equal to <n>")
    # Current index in each list in vals
    inds = [0 for _ in vals]
    greatest = []
    while len(greatest) < n:
        sub_vals = []
        # Collect values of each val at its corresponding ind in inds.
        for ind in enumerate(inds):
            if ind[1] < len(vals[ind[0]]):
                sub_vals.append(vals[ind[0]][ind[1]])
            else:
                sub_vals.append(- MaxInt() if not least else MaxInt())
        max_ind = max(range(len(sub_vals)), key=sub_vals.__getitem__)
        greatest.append(vals[max_ind][inds[max_ind]])
        inds[max_ind] += 1
    return greatest


def get_these_inds(indices: List[int], values: List[Any]) -> List[Any]:
    """Returns the values in <values> at the indices in <indices>.
    Precondition:
            For every value n in indices, 0 <= n < len(values)"""
    vals = []
    for i in indices:
        vals.append(values[i])
    return vals


def get_these_inds_matrix(indices: List[Tuple[int, int]],
                          values: List[List[Any]]) -> List[Any]:
    """Returns the values in <values> at the indices in <indices>, where
    returned_vals[0] = values[indices[0]][indices[1]]
    Precondition:
            For every value n in indices, 0 <= n < len(values)"""
    vals = []
    for tup in indices:
        vals.append(values[tup[0]][tup[1]])
    return vals


def get_cross_iteration_pattern(num_iter: int) -> List[List[Tuple[int, int]]]:
    """Returns an iteration pattern that, when performed over a square
    matrix, ensures that each column and row is iterated over the once number
    each iteration.
    Usage:
    >>> # Some square matrix, below example just shows symmetry of iteration.
    >>> matrix = [[1, 2, 3],
    ...           [4, 5, 6],
    ...           [7, 8, 9]]
    >>> for iteration in get_cross_iteration_pattern(len(matrix))
    ...     total = 0
    ...     for indices in iteration:
    ...         total += matrix[indices[0]][indices[1]]
    ...     assert total == 15
    """
    iter_p = []
    # Get first iteration, straight down the middle.
    for i in range(num_iter):
        iter_s = []
        for j in range(num_iter):
            k = j + i
            if k >= num_iter:
                k -= num_iter
            iter_s.append((j, k))
        iter_p.append(iter_s)

    return iter_p


def rev_seq(seq: Seq) -> Seq:
    """Returns a Seq with inverted directionality. If <seq> is in the 5' - 3'
    direction, then the returned seq will be in the 3' - 5' direction."""
    return Seq(str(seq)[::-1])


class MaxInt(Infinity, int):
    """An int that will always be greater than another int when compared to
    it."""
    pass


class IncompatibleSpacerError(Exception):
    """Raised when an incompatible spacer is used during some process."""
    pass
