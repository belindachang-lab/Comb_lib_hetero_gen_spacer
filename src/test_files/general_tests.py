import random
from src import hetero_spacer_generator as rsg
from src.hetero_spacer_generator.spacer_generator.random_spacer_generator import \
    gen_heterogeneity_spacers_rand, \
    get_potential_bases, \
    get_vacant_bases
import fixtures_and_helpers as fah


seqs = [
    Seq("ATCGATCG"),
    Seq("ATCG"),
    Seq("AAAAAAAA"),
    Seq("TTTTTTTT"),
    Seq(""),
    Seq("A")
]


class TestSequenceTools:
    """Test suite for sequence tools"""

    def test_get_max_complementarity(self):
        seq = seqs[1]
        assert get_max_complementarity(seq, [seqs[0]]) == 2
        assert get_max_complementarity(seq, [seqs[1]]) == 0
        assert get_max_complementarity(seq, [seqs[2]]) == 1
        assert get_max_complementarity(seq, [seqs[3]]) == 1
        assert get_max_complementarity(seq, [seqs[4]]) == 0
        assert get_max_complementarity(seq, [seqs[5]]) == 1

    def test_get_site_complementarity(self):
        assert get_site_complementarity(seqs[4], seqs[1], 1) == 0
        assert get_site_complementarity(seqs[1], seqs[0], 0) == 0
        assert get_site_complementarity(seqs[1], seqs[0], 1) == 2


def test_co_sort() -> None:
    ls = [1, 2, 3, 4, 5, 6, 7]
    lst = [2, 7, 1, 3, 5, 6, 4]
    lst1 = [2, 7, 1, 3, 5, 6, 4]
    lstc = lst[:]
    lstc1 = lst1[:]

    co_sort(lstc, lstc1)
    assert lstc == ls and lstc1 == ls


def test_co_sort_rev() -> None:
    ls = [1, 2, 3, 4, 5, 6, 7]
    lst = [2, 7, 1, 3, 5, 6, 4]
    lst1 = [2, 7, 1, 3, 5, 6, 4]
    lstc = lst[:]
    lstc1 = lst1[:]

    co_sort(lstc, lstc1, reverse=True)
    assert lstc == ls[::-1] and lstc1 == ls[::-1]


def test_co_insert() -> None:
    vals = [1, 2, 3, 5, 6]
    to_fol = [1, 1, 1, 1, 1]
    val = 4
    to_add = 0
    co_insert(vals, to_fol, val, to_add)
    assert vals == [1, 2, 3, 4, 5] and to_fol == [1, 1, 1, 0, 1]


def test_co_insert_rev() -> None:
    vals = [1, 2, 3, 5, 6][::-1]
    to_fol = [1, 1, 1, 1, 1][::-1]
    val = 4
    to_add = 0
    co_insert(vals, to_fol, val, to_add, reverse=True)
    assert vals == [2, 3, 4, 5, 6][::-1] and to_fol == [1, 1, 0, 1, 1]


def test_get_n_lowest_matrix() -> None:
    scores = [[10, 10, 10, 6],
              [10, 2, 10, 10],
              [10, 3, 10, 4],
              [5, 10, 10, 1]]
    inds = [(3, 3), (1, 1), (2, 1), (2, 3), (3, 0), (0, 3)]
    scrs = [1, 2, 3, 4, 5, 6]
    rtrn_tup = get_n_lowest_matrix(scores, 6)
    assert rtrn_tup == (inds, scrs)


def test_get_n_lowest_matrix_highest() -> None:
    scores = [[0, 0, 0, 6],
              [0, 2, 0, 0],
              [0, 3, 0, 4],
              [5, 0, 0, 1]]
    inds = [(3, 3), (1, 1), (2, 1), (2, 3), (3, 0), (0, 3)][::-1]
    scrs = [1, 2, 3, 4, 5, 6][::-1]
    assert get_n_lowest_matrix(scores, 6, highest=True) == (inds, scrs)


def test_get_lowest() -> None:
    scores = [4, 9, 1, 5, 8, 7, 2, 6, 3]
    n = 4
    inds = [2, 6, 8, 0]
    scrs = [1, 2, 3, 4]
    assert get_n_lowest(scores, n) == (inds, scrs)


def test_get_lowest_rev() -> None:
    scores = [4, 9, 1, 5, 8, 7, 2, 6, 3]
    n = 4
    inds = [1, 4, 5, 7]
    scrs = [9, 8, 7, 6]
    assert get_n_lowest(scores, n, highest=True) == (inds, scrs)


def test_get_these_inds() -> None:
    vals = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    inds = []
    for i in range(10):
        inds.append(random.choice(vals))
    new_vals = get_these_inds(inds, vals)
    for val in new_vals:
        assert val == vals[val]


def test_get_all_arrangements_basic() -> None:
    lst = get_all_arrangements(4, 4)
    assert len(lst) == 24


def test_get_cross_iteration_pattern_basic() -> None:
    matrix = [[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]]

    for iteration in get_cross_iteration_pattern(len(matrix)):
        total = 0
        for indices in iteration:
            total += matrix[indices[0]][indices[1]]
        assert (total == 15)


class TestSpacerAlignmentGen:

    def test_get_smallest_total_len_list(self):
        spacers = fah.gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        scores = get_smallest_total_len_list(spacers)

        for i in range(len(spacers)):
            assert sum(spacers[i]) == scores[i]

    def test_get_smallest_of_any_spacer_list(self) -> None:
        spacers = fah.gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        scores = get_smallest_of_any_spacer_list(spacers)

        for i in range(len(spacers)):
            assert max(spacers[i]) == scores[i]

    def test_sort_spacer_combos(self) -> None:
        spacers = fah.gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        sag = SpacerAlignmentGen(MAX_SPACER_LENGTH, NUM_HETERO)
        sag.sort_spacer_combos(spacers)
        assert sum(spacers[0]) < sum(spacers[-1]) and \
               max(spacers[0]) < max(spacers[-1])

    def test_get_all_spacer_combos_always_align(self):
        """When the heterogenity spacers are allowed to be as long as the
        heterogeneity region, there should always exist valid spacers."""
        sfm = fah.SeqFixtureManager()
        sfm.set_primers_random()
        sfm.incomplete_forward_primer.set_binding_seq("AAAAAAAAAAAA")
        sag = SpacerAlignmentGen(12, 12)
        assert sag.get_all_spacer_combos(
            sfm.incomplete_forward_primer.get_binding_seq())

    def test_always_valid_align(self) -> None:
        for i in range(fah.SMALL_SAMPLE_SIZE):
            for j in range(fah.SMALL_SAMPLE_SIZE):
                sfm = fah.SeqFixtureManager()
                sfm.for_num_hetero = i
                sfm.for_max_spacer_length = j
                sfm.set_primers_random()
                sfm.set_pot_spacers()
                for spacer in sfm.pot_forward_spacers:
                    assert fah.ensure_valid_spacers(sfm.incomplete_forward_primer,
                                                spacer, i)
                for spacer in sfm.pot_reverse_spacers:
                    assert fah.ensure_valid_spacers(sfm.incomplete_reverse_primer,
                                                spacer, i)

    def test_generates_valid_aligns(self):
        seq = Seq('GCCGGCATGGTCATGAAG')
        sag = SpacerAlignmentGen(NUM_HETERO, MAX_SPACER_LENGTH)
        primer_alignments = sag.get_all_spacer_combos(seq)
        for spacer in primer_alignments:
            assert fah.ensure_valid_spacers(seq, spacer, NUM_HETERO)

    def test_is_invalid_seq(self):
        seq = Seq('GCCGGCATGGTCATGAAG')
        sag = SpacerAlignmentGen(NUM_HETERO, MAX_SPACER_LENGTH)
        for i in range(NUM_HETERO):
            for j in range(max(1, i), NUM_HETERO):
                assert not fah.ensure_valid_spacers(seq, (0, 0, i, j), NUM_HETERO)


class TestRandomBaseSelection:
    """Test suite for methods that aid in the creation of random heterogeneity
    spacers"""

    random_tests_to_run = fah.LARGE_SAMPLE_SIZE
    hg = HeteroGen(num_hetero=3)
    hg1 = HeteroGen()
    simple_seq_arr = [['A', 'T', 'C', 'G'],
                      ['', 'A', 'T', 'C'],
                      ['', '', 'A', 'T'],
                      ['', '', '', 'A']]

    def test_gen_sequence_array(self) -> None:

        assert rsg.gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3)) \
               == self.simple_seq_arr

        for _ in range(self.random_tests_to_run):
            seq = fah.gen_random_seq(12)
            spacers = self.hg1._alignment_gen.get_all_spacer_combos(seq)
            for spacer in spacers:
                rsg.gen_sequence_array(seq, spacer)

    def test_get_vacant_bases(self) -> None:
        assert get_vacant_bases((0, 1, 2, 3)) == [
            [1, 2, 3], [2, 3],
            [3], []]

        for _ in range(self.random_tests_to_run):
            seq = fah.gen_random_seq(12)
            spacers = self.hg1._alignment_gen.get_all_spacer_combos(seq)
            for spacer in spacers:
                get_vacant_bases(spacer)

    def test_get_potential_bases(self):

        for i in range(4):
            pbs = get_potential_bases(self.simple_seq_arr, i)
            for j in range(4):
                assert self.simple_seq_arr[j][i] not in pbs

    def test_select_and_set_simple(self):

        potential_bases = ['T', 'C', 'G']
        unfilled_bases = [[1, 2, 3], []]
        column = 0
        seq_arr = self.simple_seq_arr.copy()
        for i in range(3):
            rsg.select_and_set(potential_bases, unfilled_bases,
                               column, seq_arr)
            assert len(potential_bases) == 2 - i
            assert len(unfilled_bases[0]) == 2 - i
            for row in range(1, 4):
                assert seq_arr[row][column] != 'A'

    def test_gen_heterogeneity_spacer_rand(self):
        seq = fah.gen_random_seq(12)
        spacer_combos = self.hg1._alignment_gen.get_all_spacer_combos(seq)
        if not spacer_combos:
            self.test_gen_heterogeneity_spacer_rand()
            return
        for spacer_lengths in spacer_combos:
            spacers = gen_heterogeneity_spacers_rand(seq, spacer_lengths)
            seq_arr = rsg.gen_sequence_array(seq, spacer_lengths)
            for i in range(4):
                for j in range(len(spacers[i])):
                    seq_arr[i][j] = spacers[i][j]
            assert fah.ensure_hetero_seq_arr(seq_arr, len(seq_arr[0]) - 1)


"""
    def test_remove_high_dimerisation(self) -> None:

    def test_cross_comapre(self) -> None:
    """

