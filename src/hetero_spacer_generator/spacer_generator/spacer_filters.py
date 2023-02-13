import abc
from abc import abstractmethod
from typing import Callable, List, Dict, Optional, Tuple
from Bio.Seq import Seq
from multiprocessing import Queue, Process, freeze_support
from src.hetero_spacer_generator.defaults import V, TIMING, NUM_PROCS, \
    NUM_PAIRINGS_TO_COMP
from src.hetero_spacer_generator.primer_tools import HalfSet, HeteroSeqTool, \
    MBPrimerBuilder, MaxInt, PairWiseCriterionSingle, PairwisePrimerSet, \
    PrimerSet, \
    SimpleCriterion, add, \
    eval_consecutive_complementarity, \
    eval_total_complementarity, evaluate_heterogen_binding_cross, \
    get_all_arrangements, get_cross_iteration_pattern, get_n_lowest, \
    get_n_lowest_matrix, \
    get_these_inds, remove_highest_scores, calculate_score
from src.hetero_spacer_generator.spacer_generator.criteria import EvalMBPrimer, \
    get_hetero_hetero_binding_criteria, \
    get_homo_hetero_binding_criteria
from statistics import mean

# Intervals for clearing worst performing HalfSets in PairwiseSpacerSorter
from src.hetero_spacer_generator.sequence_tools import is_degen

FIRST_INTERVAL = 2
REMOVE_HALFSET_INTERVAL = 2

# A set of four alignments, i.e. the amount by which each primer should be
# shifted to the right in a set of four to ensure heterogeneity.
SpacerAlignment = Tuple[int, int, int, int]

# A list of four sets of forward or reverse heterogeneity spacer sequences.
SpacerSet = Tuple[Seq, Seq, Seq, Seq]

if TIMING:
    import timeit


class SpacerSorter(HeteroSeqTool, abc.ABC):
    """A class that sorts spacers based on some parameters.
    _num_pairings_to_compare:
            The number of pairings of spacers to test against each other when
            evaluating pairings."""
    _num_pairings_to_compare: int

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 num_pairings_to_compare: int) -> None:
        super().__init__(max_spacer_length, num_hetero)
        self._num_pairings_to_compare = num_pairings_to_compare

    @abstractmethod
    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int) -> List[PrimerSet]:
        pass

    def set_num_pairings_to_compare(self, num_pairings: int) -> None:
        """Sets <self._num_pairings_to_compare> to <num_pairings>."""
        self._num_pairings_to_compare = num_pairings


class SortForPairwise(SpacerSorter):
    """Given sets of forward and reverse heterogeneity sequences and the
    incomplete forward and reverse primer, tries to isolate pairings of primers
    that are most likely to ensure a successful PCR.

    _for_primer:
            An incomplete forward primer.
    _rev_primer:
            An incomplete reverse primer.

    _for_seqs:
            A list of heterogeneity spacer sequences meant to pair with
            _for_primer.
    _rev_seqs:
            A list of heterogeneity spacer sequences meant to pair with
            _rev_primer.
    _initial_num_seqs:
            The initial number of _for_seqs passed.

    _single_criteria:
            A list of callables that evaluate some aspect of a spacer
            set in conjunction with some incomplete primer. Returns a list of
            scores equal to the number of spacer_sets.
    _single_criteria_weights:
            A dictionary mapping some criterion to its weight.(i.e. its relevant
            importance when evaluating primers.

    _for_single_scores:
            The scores of the forward primer sets.
    _rev_single_scores:
            The scores of the reverse primer sets.


    _pair_criteria:
            A list of callables that rate a list of forward and reverse spacers
            according to some criteria, where the first <SpacerSets> contains
            the forward spacer seqs and the first <MBPrimerBuilder> contains the
            incomplete forward primer. Returns a list of the scores.
    _pair_criteria_weights:
            A dictionary mapping some criterion to its weight.(i.e. its relevant
            importance when evaluating primers.
    _pair_criteria_weights_list:
            A list of weights where _pair_criteria_weights_list[i] is the weight
            of _pair_criteria[i]

    _pair_scores:
            A matrix mapping the indices of forward primers (i) and indices of
            reverse primers (j) to the score of that pairing
            (_pair_scores[i][j]).

    _primer_arrangements:
            A list of all potential primer arrangements, where for each
            primer_arrangement in _primer_arrangements and for any set of
            forward and reverse primers forward_primers[i] pairs with
            reverse_primers[primer_arrangement[i]].

    _for_halfsets:
            HalfSets constructed from the best forward SpacerSets.
    _rev_halfsets:
            HalfSets constructed from the best reverse SpacerSets.
    _fullsets:
            Full sets of primers constructed from the forward and reverse
            halfsets.

    _num_procs:
            The number of processes to spawn when multithreading.
        """

    _degen: bool

    _primer_evaluator: EvalMBPrimer

    _for_primer: MBPrimerBuilder
    _rev_primer: MBPrimerBuilder

    _for_seqs: List[SpacerSet]
    _rev_seqs: List[SpacerSet]
    _initial_num_seqs: int

    _single_criteria: List[PairWiseCriterionSingle]
    _single_criteria_weights: Dict[PairWiseCriterionSingle, int]

    _for_single_scores: List[int]
    _rev_single_scores: List[int]

    _pair_criteria: List[SimpleCriterion]
    _pair_criteria_weights: Dict[SimpleCriterion, int]
    _pair_criteria_weights_list: List[int]

    _for_halfsets: List[HalfSet]
    _rev_halfsets: List[HalfSet]

    _fullsets: List[List[PairwisePrimerSet]]

    _num_procs: int

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 num_pairings_to_comp: int = NUM_PAIRINGS_TO_COMP,
                 degen: bool = None) -> None:
        """May initialise a functionally empty class. Build method needs to be
        called in order to ensure proper method functionality."""
        super().__init__(max_spacer_length, num_hetero, num_pairings_to_comp)
        self._degen = degen
        self._primer_evaluator = EvalMBPrimer(max_spacer_length, num_hetero,
                                              self._degen)
        self._for_primer = MBPrimerBuilder()
        self._rev_primer = MBPrimerBuilder()
        self._for_seqs = []
        self._rev_seqs = []
        self._for_single_scores = []
        self._rev_single_scores = []
        self._pair_criteria = []
        self._num_procs = NUM_PROCS
        self._build_criteria()

    def _do_degen_check(self, primer: MBPrimerBuilder) -> None:
        """Checks the sequences to be used in this filter for degeneracy,
        adjusts class attributes accordingly."""
        for seq in primer:
            if is_degen(seq):
                self._degen = True
                return
        self._degen = False
        return

    def _build_full(self, max_spacer_length: int, num_hetero: int,
                    for_primer: MBPrimerBuilder, rev_primer: MBPrimerBuilder,
                    for_seqs: List[SpacerSet], rev_seqs: List[SpacerSet]
                    ) -> None:
        """Initialises the classes attributes to the given values.
        Precondition:
                len(for_seqs) == len(rev_seqs)"""
        super().__init__(max_spacer_length, num_hetero,
                         self._num_pairings_to_compare)
        if not self._degen:
            self._do_degen_check(for_primer)
            self._do_degen_check(rev_primer)
        self._primer_arrangements = get_all_arrangements(4, 4)
        self._for_primer = for_primer
        self._rev_primer = rev_primer
        self._for_seqs = for_seqs
        self._rev_seqs = rev_seqs
        self._initial_num_seqs = len(for_seqs)
        self._for_single_scores = []
        self._rev_single_scores = []
        self._for_halfsets = []
        self._rev_halfsets = []
        for i in range(len(for_seqs)):
            self._rev_single_scores.append(1)
            self._for_single_scores.append(1)

    # INSERT NEW CRITERIA HERE
    def _build_criteria(self) -> None:
        """Constructs the criteria attributes with some Callables"""
        self._single_criteria, self._single_criteria_weights = \
            get_homo_hetero_binding_criteria(self._primer_evaluator)

        self._pair_criteria, self._pair_criteria_weights = \
            get_hetero_hetero_binding_criteria(self._primer_evaluator)
        self._pair_criteria_weights_list = []
        for criterion in self._pair_criteria:
            self._pair_criteria_weights_list.append(
                self._pair_criteria_weights[criterion])

    # Begin methods for finding seqs with lowest homo-dimerisation potential
    # from a pre-existing set of primers.

    def _evaluate_scores_single(self):
        """Evaluates the scores of both the forward and reverse primers for each
         criterion in <self._single_criteria>"""
        for criterion in self._single_criteria:
            self._get_scores_single(self._for_primer, self._for_seqs, criterion,
                                    self._for_single_scores)
            self._get_scores_single(self._rev_primer, self._rev_seqs, criterion,
                                    self._rev_single_scores)

    def _get_scores_single(self, primer: MBPrimerBuilder,
                           spacer_seqs: List[SpacerSet],
                           criterion: PairWiseCriterionSingle,
                           scores: List[int]) -> None:
        """Updates the score of each SpacerSet in <scores> using the <criterion>
        method."""
        weight = self._single_criteria_weights[criterion]
        for i in range(len(spacer_seqs)):
            set_scores = []
            for seq in spacer_seqs[i]:
                primer.set_heterogen_seq(seq)
                comp_primer = primer.get_mbprimer()
                set_scores.append(criterion(comp_primer))
            scores[i] += calculate_score(set_scores) * weight

    def _sort_and_trim(self):
        """Removes worst scoring forward and reverse primers and their scores
        from their respective attributes."""
        best_for_inds, best_for_scores = \
            get_n_lowest(self._for_single_scores,
                         self._num_pairings_to_compare)
        best_rev_inds, best_rev_scores = \
            get_n_lowest(self._rev_single_scores,
                         self._num_pairings_to_compare)

        self._for_seqs, self._for_single_scores = \
            get_these_inds(best_for_inds, self._for_seqs), best_for_scores
        self._rev_seqs, self._rev_single_scores = \
            get_these_inds(best_rev_inds, self._rev_seqs), best_rev_scores

    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int, degen: bool = None,
                                    num_procs: int = NUM_PROCS) \
            -> List[PairwisePrimerSet]:
        """Returns the PrimerSets best suited for pairwise PCR for
        metabarcoding."""
        self._build_full(self._max_spacer_length, self._num_hetero,
                         incomplete_forward_primer, incomplete_reverse_primer,
                         forward_spacers, reverse_spacers)
        self._num_procs = num_procs
        if TIMING:
            single_scoring = self._evaluate_scores_single
            double_scoring = self._score_primer_sets

            self._primer_evaluator.start_counting_comparisons()
            sst = timeit.timeit(single_scoring, number=1)
            count = self._primer_evaluator.stop_counting_comparisons()
            print("Primer sets evaluated individually (2x{num_primers}) in "
                  "{time:.2f} seconds using "
                  "{cnt} sequence comparisons.".
                  format(time=sst, cnt=count, num_primers=len(self._for_seqs)))

            self._sort_and_trim()

            self._primer_evaluator.start_counting_comparisons()
            dst = timeit.timeit(double_scoring, number=1)
            count = self._primer_evaluator.stop_counting_comparisons()
            print("Primer evaluated together ({len1}x{len2}) in {time:.2f} "
                  "seconds using {cnt} sequence comparisons.".
                  format(time=dst, cnt=count, len1=len(self._for_halfsets),
                         len2=len(self._rev_halfsets)))

        else:
            self._evaluate_scores_single()
            self._sort_and_trim()
            self._score_primer_sets()

        if self._num_pairings_to_compare == 1:
            # Randomise the pairing of the best set. Produces entirely
            # un-optimised result.
            self._fullsets[0][0].randomise_optimal_pairing()
        self._print_or_csv(V, False)

        return self._get_lowest_scoring_sets(num_to_return)

    @staticmethod
    def get_csv_formatting_str() -> str:
        """Returns a string with formatting of csv returned by
        self._print_or_csv."""
        return 'sets evaluated for homogeneity (x2), ' \
               'pairs evaluated for heterogeneity(~^2), ' \
               'filtered forward homo avg,' \
               'filtered reverse homo avg,' \
               'forward homo min,' \
               'reverse homo min,' \
               'hetero avg, ' \
               'hetero min, ' \
               'best set forward homogeneity, ' \
               'best set reverse homogeneity, ' \
               'best set heterogeneity'

    def _print_attributes(self) -> None:
        """Prints the attributes of the current *completed* filter results."""
        self._print_or_csv(True, False)

    def get_csv(self) -> Tuple[float]:
        """Returns the attributes of the current *completed* filter results in
        csv format."""
        return self._print_or_csv(False, True)

    def _print_or_csv(self, prnt: bool, csv: bool) -> Optional[Tuple[float]]:
        """Iff csv, returns a str rep of a tuple of various score attributes.
        Iff <prnt>, prints a string representation of the results."""
        if prnt or csv:
            # Compile the values required for csv string or print.
            sets_eval_homo = self._initial_num_seqs
            pairs_eval_hetero = self._num_pairings_to_compare
            forward_homo_avg = mean(self._for_single_scores)
            reverse_homo_avg = mean(self._rev_single_scores)
            forward_homo_min = min(self._for_single_scores)
            reverse_homo_min = min(self._rev_single_scores)

            hetero_sum = 0
            hetero_num = 0
            hetero_min = MaxInt()
            min_set_f = -1
            min_set_r = -1
            # For each FullSet if it has been scored, include it in the average.
            for f in range(self._num_pairings_to_compare):
                for r in range(self._num_pairings_to_compare):
                    fs = self._fullsets[f][r]
                    # Select fullset, update values if neccesaru
                    if fs.has_been_scored():
                        hetero_sum += fs.get_score()
                        hetero_num += 1
                        if fs.get_score() < hetero_min:
                            hetero_min = fs.get_score()
                            min_set_f = f
                            min_set_r = r

            hetero_avg = hetero_sum / hetero_num

            best_for_homo = self._for_single_scores[min_set_f]
            best_rev_homo = self._rev_single_scores[min_set_r]
            best_hetero = hetero_min

        # Print out the string describing the findings from this run.
            if prnt:
                # Plurality string
                msg = 'Hetero Comp: {cc} - Homo Comp: {hc}\n'.format(
                    cc=sets_eval_homo, hc=pairs_eval_hetero
                )
                msg += "For Homo Avg: {fa:.2f} - Rev Homo Avg: {ra:.2f}\n".format(
                    fa=forward_homo_avg,
                    ra=reverse_homo_avg
                )
                msg += "For Homo Min: {fa} - Rev Homo Min: {ra}\n".format(
                    fa=forward_homo_min,
                    ra=reverse_homo_min
                )


                msg += "Hetero Avg: {ha:.3f}\n".format(
                    ha=hetero_avg
                )
                msg += "Best primer set:\n"

                msg += "Hetero: {hs} - For Homo: {fh} - Rev Homo {rh}\n" \
                        .format(hs=best_hetero, fh=best_for_homo,
                                rh=best_rev_homo
                                )
                print(msg)

            # If a CSV was requested, return it.
            if csv:
                tup = (
                    sets_eval_homo,
                    pairs_eval_hetero,
                    forward_homo_avg,
                    reverse_homo_avg,
                    forward_homo_min,
                    reverse_homo_min,
                    hetero_avg,
                    hetero_min,
                    best_for_homo,
                    best_rev_homo,
                    best_hetero
                )
                return tup
        return None

    def _construct_half_sets(self) -> None:
        """Returns (forward_halfset, reverse_halfset) constructed from the
        forward and reverse primers and the sets of heterogeneity spacer
        sequences."""
        for for_primer_seqs in self._for_seqs:
            self._for_halfsets.append(HalfSet(self._for_primer,
                                              for_primer_seqs))
        for rev_primer_seqs in self._rev_seqs:
            self._rev_halfsets.append(HalfSet(self._rev_primer,
                                              rev_primer_seqs))
        return

    def _construct_full_set_matrix(self) -> None:
        """Constructs a matrix of PairwisePrimerSets from the lists of forward
        and reverse HalfSets."""
        self._fullsets = []
        for f in range(self._num_pairings_to_compare):
            self._fullsets.append([])
            for r in range(self._num_pairings_to_compare):
                self._fullsets[f].append(PairwisePrimerSet(
                    self._for_halfsets[f], self._rev_halfsets[r]))
        return

    def _get_fullset_matrix_string(self) -> str:
        """Returns a string representation of the current fullset matrix and the
        scores it contains."""

        tab_len = 5
        def app_get_ts(s: str, nt: int = 1) -> str:
            """Appends tabs as neccesary so that the string fits into a 2 tab/
            column matrix string."""

            ns = nt * tab_len - len(s)
            return s + ' ' * ns

        rtrn_str = app_get_ts('', nt=3)
        # Forward
        for f, hs in enumerate(self._for_halfsets):
            rtrn_str += app_get_ts("FS{ind}".format(ind=f))
        rtrn_str += '\n'
        rtrn_str += app_get_ts('', nt=3)
        for f, hs in enumerate(self._for_halfsets):
            rtrn_str += app_get_ts(str('{score:.1f}'.format(score=hs.get_avg())))
        rtrn_str += '\n'
        for r, rhs in enumerate(self._rev_halfsets):
            rtrn_str += app_get_ts("RS{ind}".format(ind=r), nt=1)
            rtrn_str += app_get_ts(' {score:.1f}'.format(
                score=rhs.get_avg()), nt=2)
            for f, fhs in enumerate(self._for_halfsets):
                # It's been scored, add its score.
                if self._fullsets[f][r].has_been_scored():
                    rtrn_str += app_get_ts('{score}'.format(
                        score=self._fullsets[f][r].get_score()))
                # Unscored and row removed. Pass.
                elif not (fhs.is_active() and rhs.is_active()):
                    rtrn_str += app_get_ts('%')
                # Not yet scored, but still included.
                else:
                    rtrn_str += app_get_ts('-')
            rtrn_str += '\n'
        return rtrn_str

    def _get_row_average(self, row: int) -> float:
        """Gets the average of the scored fullsets in some <row> of
        <self._fullsets>."""
        sum_fullsets = 0
        num = 0
        for i in range(self._num_pairings_to_compare):
            if self._fullsets[row][i].has_been_scored():
                sum_fullsets += self._fullsets[row][i].get_min_pairing_score()
                num += 1
        return sum_fullsets / num

    def _get_col_average(self, col: int) -> float:
        """Gets the average of the scored fullsets in some column (<col>) of
        <self._fullsets>."""
        sum_fullsets = 0
        num = 0
        for i in range(self._num_pairings_to_compare):
            if self._fullsets[i][col].has_been_scored():
                full_set = self._fullsets[i][col]
                sum_fullsets += full_set.get_min_pairing_score()
                num += 1
        return sum_fullsets / num

    def _update_averages(self) -> None:
        """Updates the averages of the HalfSets."""
        for f in range(self._num_pairings_to_compare):
            if self._for_halfsets[f].is_active():
                row_avg = self._get_row_average(f)
                self._for_halfsets[f].set_avg(row_avg)

        for r in range(self._num_pairings_to_compare):
            if self._for_halfsets[r].is_active():
                col_avg = self._get_col_average(r)
                self._rev_halfsets[r].set_avg(col_avg)
        return

    def _remove_worst_performing(self) -> None:
        """Removes the worst performing HalfSet from both the forward and
        reverse halfsets. """
        for_scores = []
        rev_scores = []
        for i in range(self._num_pairings_to_compare):
            if self._for_halfsets[i].is_active():
                for_scores.append(self._for_halfsets[i].get_avg())
            else:
                for_scores.append(0)
            if self._rev_halfsets[i].is_active():
                rev_scores.append(self._rev_halfsets[i].get_avg())
            else:
                rev_scores.append(0)
        worst_for_index = get_n_lowest(for_scores, 1, highest=True)[0][0]
        worst_rev_index = get_n_lowest(rev_scores, 1, highest=True)[0][0]
        self._for_halfsets[worst_for_index].deactivate()
        self._rev_halfsets[worst_rev_index].deactivate()
        return

    def _score_primer_sets(self) -> None:
        """Creates and scores all combinations of
        <self._num_pairings_to_compare> the forward and reverse primer sets."""
        self._construct_half_sets()
        self._construct_full_set_matrix()
        # Symmetric iteration pattern.
        iter_pttrn = get_cross_iteration_pattern(self._num_pairings_to_compare)
        if self._num_procs > 1:
            # For communicating indices of FullSets to be evaluated to child
            # processes and receiving completed FullSets.
            # format: Tuple(int, int)
            inds = Queue()
            # format: Tuple(int, int, FullSet)
            completed = Queue()
            processes = []
            for _ in range(self._num_procs):
                processes.append(Process(target=self._proc_eval_fullset_pair,
                                         args=(inds, completed)))
            for process in processes:
                process.start()

        # Iterate through the full set. Is costly to evaluate set with criteria,
        # so continually remove the worst performing ones.
        for i in range(self._num_pairings_to_compare):

            # Remove worst scoring set after a given interval
            if i - FIRST_INTERVAL >= 0 \
                    and (i - FIRST_INTERVAL) % REMOVE_HALFSET_INTERVAL == 0:
                self._update_averages()
                self._remove_worst_performing()

            to_collect = 0
            # Delegate responsibility to the processes
            for j in range(self._num_pairings_to_compare):
                f, r = iter_pttrn[i][j]
                # If sets are still active, do evaluation.
                if self._for_halfsets[f].is_active() \
                        and self._rev_halfsets[r].is_active():
                    if self._num_procs > 1:
                        to_collect += 1
                        inds.put((f, r))
                    else:
                        self._fullsets[f][r].apply_criteria(
                            self._pair_criteria,
                            self._pair_criteria_weights_list)

            if self._num_procs > 1:
                # Collect FullSets
                for _ in range(to_collect):
                    f, r, fullset = completed.get()
                    self._fullsets[f][r] = fullset

            if V:
                percent_complete = (i + 1) / self._num_pairings_to_compare * 100
                self._update_averages()
                print("{prcnt:.2f}% complete...".
                      format(prcnt=percent_complete))
                print(self._get_fullset_matrix_string())

        if self._num_procs > 1:
            # Add stop commands to queue.
            for _ in range(self._num_procs):
                inds.put((-1, -1))

    def _proc_eval_fullset_pair(self, q_in: Queue, q_out: Queue) -> None:
        """Evaluates the fullset at the given position."""
        if __name__ == '__main__':
            freeze_support()
        f, r = q_in.get()
        while f >= 0:
            self._fullsets[f][r].apply_criteria(
                self._pair_criteria,
                self._pair_criteria_weights_list)
            q_out.put((f, r, self._fullsets[f][r]))
            f, r = q_in.get()

    def _get_lowest_scoring_sets(self, num_sets: int) \
            -> List[PairwisePrimerSet]:
        """Constructs a matrix containing the scores of the now evaluated
        PairwisePrimerSets."""
        scores = []
        for f in range(self._num_pairings_to_compare):
            scores.append([])
            for r in range(self._num_pairings_to_compare):
                if self._fullsets[f][r].has_been_scored():
                    scores[f].append(self._fullsets[f][r].get_min_pairing_score())
                else:
                    scores[f].append(MaxInt())
        lowest_scoring_inds = get_n_lowest_matrix(scores, num_sets)[0]
        lowest_scoring = []
        for ind in lowest_scoring_inds:
            f, r = ind
            lowest_scoring.append(self._fullsets[f][r])

        return lowest_scoring


def remove_high_consec_complementarity(spacers: List[SpacerSet],
                                       incomplete_primer: MBPrimerBuilder,
                                       num_to_keep: int) -> None:
    """Removes the spacers from <spacers> that have the greatest number of
    consecutive complementary bases."""
    # Map consecutive complementarity of spacers to their index in
    # <spacers>.
    complementarity_to_index = {}
    for i in range(len(spacers)):
        complement = eval_consecutive_complementarity(incomplete_primer,
                                                      spacers[i])
        add(complementarity_to_index, complement, i)

    # Remove the spacers with highest complementarity.
    remove_highest_scores(spacers, complementarity_to_index, num_to_keep)


def remove_high_dimer_complementarity(spacers: List[SpacerSet],
                                      incomplete_primer: MBPrimerBuilder,
                                      num_to_keep: int) -> None:
    """Removes the spacers from <spacers> with the highest self
    complementarity with <incomplete_primer>. Retains the <num_to_keep>
    spacers with the lowest self complementarity"""
    # Map complementarity of spacers to their index in <spacers>.
    complementarity_to_index = {}
    for i in range(len(spacers)):
        complement = eval_total_complementarity(incomplete_primer,
                                                spacers[i])
        add(complementarity_to_index, complement, i)

    # Remove the spacers with highest complementarity.
    remove_highest_scores(spacers, complementarity_to_index, num_to_keep)


def cross_compare(incomplete_forward_primer: MBPrimerBuilder,
                  incomplete_reverse_primer: MBPrimerBuilder,
                  forward_spacer_seqs: List[SpacerSet],
                  reverse_spacer_seqs: List[SpacerSet],
                  num_to_return: int) -> List[PrimerSet]:
    """Compares various combinations of spacer sets, raking them from lowest
    inter-set heterogeneity spacer / primer complementarity"""
    forward_sets = []
    reverse_sets = []
    # Convert sets of spacers & incomplete primers into sets of complete
    # primers. Assumes that len(forward_spacer_seqs) == len(
    # reverse_spacer_seqs)
    for i in range(len(forward_spacer_seqs)):
        forward_sets.append(HalfSet(incomplete_forward_primer,
                                    forward_spacer_seqs[i]))
        reverse_sets.append(HalfSet(incomplete_reverse_primer,
                                    reverse_spacer_seqs[i]))

    scores = []
    # Create a matrix containing the scores between any 2 of the forward and
    # reverse primers.
    for f in range(len(forward_sets)):
        scores.append([])
        for r in range(len(reverse_sets)):
            scores[f].append(
                evaluate_heterogen_binding_cross(forward_sets[f],
                                                 reverse_sets[r]))

    best_sets = []
    min_combos = get_n_lowest_matrix(scores, num_to_return)[0]
    for i in range(num_to_return):
        p_set = PrimerSet(forward_sets[min_combos[i][0]],
                          reverse_sets[min_combos[i][1]])
        best_sets.append(p_set)
    return best_sets


class SortForSimultaneous(SpacerSorter):
    """Sorts spacers according to some parameters. Returns sets of primers that
    would bind minimally when all used simultaneously in PCR.

    === Private Attributes ===
    _number_to_cross_compare:
            The number of forward and reverse spacers to compare when evaluating
            primer hetero-dimer formation.
    _criteria:
            A list of functions or methods that sort a list of potential spacer
             sequences <SpacerSets> for some <MBPrimerBuilder> according
             to some metric, and leave some num <int> of the best scoring sets
             of spacer sequences in the list."""

    _criteria: List[Callable[[List[SpacerSet], MBPrimerBuilder, int], None]]

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 rigour: int = 1) -> None:
        """Initialises the attributes to the values specified. Sets the sampling
        size attributes (self._random_per_align & self._number_to_cross_compare)
        in accordance with the rigour.

        Note: increases to <rigour> will exponentially increase the runtime of
        self._get_hetero_seqs, but may increase the quality of the heterogeneity
        sequences produced."""
        super().__init__(max_spacer_length, num_hetero,
                         NUM_PAIRINGS_TO_COMP * rigour)
        self._criteria = [
            remove_high_dimer_complementarity,
            remove_high_consec_complementarity
        ]

    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences.
        Note:
            Will mutate forward and reverse spacers."""

        # Filter out all but prcnt_to_keep% potential heterogeneity seqs.
        prcnt_to_keep = (self._num_pairings_to_compare /
                         len(forward_spacers) * 100)
        self._filter_spacer_sets(incomplete_forward_primer,
                                 incomplete_reverse_primer,
                                 forward_spacers, reverse_spacers,
                                 prcnt_to_keep)

        # Find the best combinations of forward and reverse primers with the
        # least spacer region complementarity.
        return cross_compare(incomplete_forward_primer,
                             incomplete_reverse_primer,
                             forward_spacers,
                             reverse_spacers, num_to_return)

    def _filter_spacer_sets(self, incomplete_forward_primer: MBPrimerBuilder,
                            incomplete_reverse_primer: MBPrimerBuilder,
                            forward_spacer_seqs: List[SpacerSet],
                            reverse_spacer_seqs: List[SpacerSet],
                            percent_to_keep: float) -> None:
        """Removes all but <percent_to_keep>% of elements in
        <forward_spacer_seqs> and <reverse_spacer_seqs>. Removes elements based
        on a set of criteria methods in <self.criteria>.
        Note:
            Rounds percentage down, so may return less spacers than expected."""
        # Calculate the amounts by which to decrement the number of spacers each
        # time
        org_len = len(forward_spacer_seqs)
        percent_to_remove_every_iter = ((100 - percent_to_keep) /
                                        len(self._criteria))
        num_to_remove_each_iter = org_len * percent_to_remove_every_iter / 100
        num_to_keep = org_len
        for criterion in self._criteria:
            num_to_keep -= num_to_remove_each_iter
            criterion(forward_spacer_seqs, incomplete_forward_primer,
                      int(num_to_keep))
            criterion(reverse_spacer_seqs, incomplete_reverse_primer,
                      int(num_to_keep))
