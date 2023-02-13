import pytest
from seq_alignment_analyser.best_primers import *
from test_files.fixtures_and_helpers import all_unique

from test_files.fixtures_and_helpers import get_msa

seqn = ("GGATTCCCCGTCAACTTCCTCACGCTGTACGTCACAATCGAACACAAGAA"
        "GCTACGCTCGCCTCTCAACTACATCCTGCTCAATCTGGCCGTGGCCGACC"
        "TCTTCATGGTGATCGGCGGCTTCACTACCACGATGTGGACCTCGCTCAAC"
        "GGCTACTTCGTCTTCGGCCGCGTGGGCTGCAACATCGAGGGCTTCTTCGC"
        "CACCCTGGGCGGTGAGATCGCGCTCTGGTCCCTGGTGGTGCTGTCCATGG"
        "AGAGGTGGATCGTCGTCTGCAAGCCCATGAGCAACTTCCGCTTCGGTGAG"
        "AACCACGCCGTTATGGGCGTCGCATTCTCATGGGTGATGGCCTCTGCTTG"
        "TGCCGTGCCTCCCCTGGTCGGCTGGTCCCGTTACATCCCTGAGGGCATGC"
        "AGTGCTCATGCGGAATTGACTACTACACACGCGCCGAAGGCTTCAACAAC"
        "GAGTCCTTTGTCATCTACATGTTCATTGTCCACTTCACCATCCCCCTGAT"
        "CATCATCACCTTCTGCTACGGCCGCCTGGTCTGCACCGTCAAGGAGGCCG"
        "CCGCCCAGCAGCAGGAGTCCGAGACCACCCAAAGGGCTGAGCGTGAGGTC"
        "ACCCGCATGGTCATCATCATGTTCGTTGCCTTCTTGGTATGCTGGGTGCC"
        "CTATGCAAGTGTAGCTTGGTATATCTTCACGCACCAGGGCAGTGAGTTCG"
        "GTCCGGTCTTCATGACCATTCCG")


@pytest.fixture
def consensus() -> str:
    """A generic DNA sequence, extracted from a rhodopsin ORF."""
    return seqn


# len = 723
# max_ind = 722

@pytest.mark.parametrize('gen_al_5p,gen_al_len', [
    # ([], []),
    # ([], []),
    ([717], [6]),
    (list(range(0, 700)), [20]),
    ([40, 41], [30, 20])
])
class TestHomoSeqIterator:
    """Tests that HomoSeqIterator produces expected results"""

    def test_all_unique(self, consensus, gen_al_5p, gen_al_len) -> None:
        """Tests that every sequence returned is unique. Assumes no repeats in
        <gen_consensus> of lengths <gen_al_len>."""
        iter = HomoSeqIterator(consensus, gen_al_5p, gen_al_len)
        rslts = [it for it in iter]
        assert all_unique(rslts)

    def test_corr_len(self, consensus, gen_al_5p, gen_al_len) -> None:
        """Tests that the iterator produces the expected number of sequences."""
        iter = HomoSeqIterator(consensus, gen_al_5p, gen_al_len)
        rslts = [it for it in iter]
        assert len(rslts) == len(gen_al_5p) * len(gen_al_len)

    def test_correct_sequences(self, consensus, gen_al_5p,
                               gen_al_len) -> None:
        """Tests that the correct sequences are returned."""
        iter = HomoSeqIterator(consensus, gen_al_5p, gen_al_len)
        rslts = [it for it in iter]
        for s5p in gen_al_5p:
            for length in gen_al_len:
                assert consensus[s5p:s5p + length] in rslts

    def test_termination(self, consensus, gen_al_5p, gen_al_len) -> None:
        """Tests that the iterator terminates as expected"""
        iter = HomoSeqIterator(consensus, gen_al_5p, gen_al_len)
        num_iters = len(gen_al_5p) * len(gen_al_len)
        for _ in range(num_iters):
            assert iter.__next__() is not None
        with pytest.raises(StopIteration):
            iter.__next__()


all_len = list(range(40, len(seqn) + 1))


@pytest.mark.parametrize('f_al_5p,r_al_5p,f_al_len,r_al_len,a_al_len',
                         [
                             # ([],[],[],[],[]),
                             ([0], [722], [20], [20], [723]),
                             ([100], [200], [10], [10], [101]),
                             ([0, 1, 2, 3], [722, 721, 720, 719], [20, 21, 22],
                              [41, 50, 60, 100], all_len),
                             ([0, 1, 2, 3], [722, 721, 720, 719], [20], [20],
                              all_len)
                         ], scope='class')
class TestHeteroSeqIteratorCompletion:
    """Tests that HomoSeqIterator produces expected results when working with
    primer sets that do not have amplicon size restrictions."""

    def test_expected_num_results(self, consensus, f_al_5p, r_al_5p, f_al_len,
                                  r_al_len, a_al_len) -> None:
        """Test that the iterator returns the appropriate number of results."""
        hsi = HeteroSeqIterator(consensus, f_al_5p, r_al_5p, f_al_len, r_al_len,
                                a_al_len)
        rslts = [rslt for rslt in hsi]
        num_f_seqs = len(f_al_5p) * len(f_al_len)
        num_r_seqs = len(r_al_5p) * len(r_al_len)
        assert len(rslts) == num_f_seqs * num_r_seqs

    def test_unique_results(self, consensus, f_al_5p, r_al_5p, f_al_len,
                            r_al_len, a_al_len) -> None:
        """Asserts that all results are unique."""
        hsi = HeteroSeqIterator(consensus, f_al_5p, r_al_5p, f_al_len, r_al_len,
                                a_al_len)
        rslts = [rslt for rslt in hsi]
        assert all_unique(rslts)

    def test_expected_seqs(self, consensus, f_al_5p, r_al_5p, f_al_len,
                           r_al_len, a_al_len) -> None:
        """Tests that the sequences are selected as expected."""
        hsi = HeteroSeqIterator(consensus, f_al_5p, r_al_5p, f_al_len, r_al_len,
                                a_al_len)
        rslts = [rslt for rslt in hsi]
        for f5p in f_al_5p:
            for r5p in r_al_5p:
                for flen in f_al_len:
                    for rlen in r_al_len:
                        # + 1 for the reverse since the primer starts at the
                        # binding region.
                        rev_bind = Seq(consensus[r5p - rlen + 1: r5p + 1])
                        rev_bind = rev_bind.reverse_complement()
                        rev_bind = str(rev_bind)
                        tup = (consensus[f5p:f5p + flen],
                               rev_bind)
                        assert tup in rslts


@pytest.mark.parametrize('f_al_5p,r_al_5p,f_al_len,r_al_len,a_al_len',
                         [
                             (list(range(99, 150)), list(range(399, 450)),
                              [20, 40, 30], [25, 30, 35], list(range(285, 316)))
                             # ([],[],[],[],[])
                             # ([],[],[],[],[])
                             # ([],[],[],[],[])
                             # ([],[],[],[],[])
                             # ([],[],[],[],[])
                         ], scope='class')
class TestHeteroSeqIteratorNoCompletion:
    """Tests for hetero sequences where not all combinations of binding
    sequences are valid."""

    def test_correct_num(self, consensus, f_al_5p, r_al_5p, f_al_len,
                         r_al_len, a_al_len) -> None:
        """Tests that the correct number of results are produced"""
        hsi = HeteroSeqIterator(consensus, f_al_5p, r_al_5p, f_al_len, r_al_len,
                                a_al_len)
        rslts = [rslt for rslt in hsi]
        num_pos = 0
        for f5p in f_al_5p:
            for r5p in r_al_5p:
                for alen in a_al_len:
                    if f5p + alen - 1 == r5p:
                        num_pos += len(f_al_len) * len(r_al_len)
                        break

        assert len(rslts) == num_pos


def test_find_average_conservation() -> None:
    """Tests that a binding seq scorer is able to properly evaluate the average
    conservation in an alignment."""
    msa = get_msa('known_conservation')
    f_hss = HomoSeqScorer(msa, [1], [1], ['a'], 4, False)
    r_hss = HomoSeqScorer(msa, [1], [1], ['a'], 4, True)

    # First and last four bases respectively.
    f1 = BindingParams(0, 4)
    r1 = BindingParams(0, 4)
    r2 = BindingParams(15, 4)

    f_hss._find_average_conservation(f1, False)
    f_hss._find_average_conservation(r2, False, dir='r')
    r_hss._find_average_conservation(r1, True)

    assert f1.get_mean_conservation() == 25
    assert r1.get_mean_conservation() == 100
    assert r2.get_mean_conservation() == 100

def weight_formula(c, t, s1, s2) -> int:
    return (c * CONSEC_WEIGHT + t * TOTAL_WEIGHT) / (len(s1) + len(s2))

def test_get_dimer_score_simple() -> None:
    """Tests that the program is succesfully able to locate dimer structures."""
    seq1 = 'AAAA'
    seq2 = 'AAAA'
    assert get_dimer_score(seq1, seq2) == weight_formula(0, 0, seq1, seq2)

def test_get_dimer_score_basic_2() -> None:
    """Tests that the program is succesfully able to locate dimer structures."""
    s1 = 'TTAA'
    s2 = 'TTAA'
    assert get_dimer_score(s1, s2) == weight_formula(4, 4, s1, s2)


def test_get_dimer_score_basic_3() -> None:
    """Tests that the program is succesfully able to locate dimer structures."""
    s1 = 'ATATATATTTTTTTTTTTTTTTTTGGG'
    s2 = 'CCCTTTTTCTCTCT'
    assert get_dimer_score(s1, s2) == weight_formula(3, 4, s1, s2)

def test_basic_operation() -> None:
    """Tests that the program is successfully able to locate dimer structures."""
    msa = get_msa('known_conservation')
    amp_lens = [10]

    f_5p = [3]
    r_5p = [12]

    f_len = [3]
    r_len = [3]

    f_a = ['AAAA']
    r_a = ['AAAA']

    # Basic HomoSeqScorer functionality.
    f_hos = HomoSeqScorer(msa, f_5p, f_len, f_a, 5, False)
    r_hos = HomoSeqScorer(msa, r_5p, r_len, r_a, 5, True)

    f_hos.find_best()
    r_hos.find_best()

    assert len(f_hos.get_best_params()) == len(r_hos.get_best_params()) == 1

    hes = HeteroSeqScorer(msa, f_5p, r_5p, f_len, r_len, amp_lens, f_a, r_a, 5)

    hes.find_best()

    assert len(hes.get_best_params()) == 1

def test_basic_operation_2() -> None:
    """Tests that the program is successfully able to locate dimer structures."""
    msa = get_msa('known_conservation')
    amp_lens = [10]

    f_5p = [1, 2, 3]
    r_5p = [10, 11, 12]

    f_len = [3]
    r_len = [3]

    f_a = ['AAAA']
    r_a = ['AAAA']

    # Basic HomoSeqScorer functionality.
    f_hos = HomoSeqScorer(msa, f_5p, f_len, f_a, 5, False)
    r_hos = HomoSeqScorer(msa, r_5p, r_len, r_a, 5, True)

    f_hos.find_best()
    r_hos.find_best()

    assert len(f_hos.get_best_params()) == len(r_hos.get_best_params()) == 3

    hes = HeteroSeqScorer(msa, f_5p, r_5p, f_len, r_len, amp_lens, f_a, r_a, 5)

    hes.find_best()

    assert len(hes.get_best_params()) == 3


def test_basic_operation_3() -> None:
    """Tests that the program is successfully able to locate dimer structures."""
    msa = get_msa('known_conservation')
    amp_lens = [8, 9, 10, 11, 12]

    f_5p = [1, 2, 3]
    r_5p = [10, 11, 12]

    f_len = [3]
    r_len = [3]

    f_a = ['AAAA']
    r_a = ['AAAA']

    # Basic HomoSeqScorer functionality.
    f_hos = HomoSeqScorer(msa, f_5p, f_len, f_a, 5, False)
    r_hos = HomoSeqScorer(msa, r_5p, r_len, r_a, 5, True)

    f_hos.find_best()
    r_hos.find_best()

    assert len(f_hos.get_best_params()) == len(r_hos.get_best_params()) == 3

    hes = HeteroSeqScorer(msa, f_5p, r_5p, f_len, r_len, amp_lens, f_a, r_a, 11)

    hes.find_best()

    assert len(hes.get_best_params()) == 9


def test_basic_operation_4() -> None:
    """Tests that the program is successfully able to locate dimer structures."""
    msa = get_msa('known_conservation')
    amp_lens = [8, 9, 10, 11, 12, 13, 14, 15, 16]

    f_5p = [0, 1, 2]
    r_5p = [13, 14, 15]

    f_len = [1, 2, 3]
    r_len = [1, 2, 3]

    f_a = ['AAAA']
    r_a = ['AAAA']

    # Basic HomoSeqScorer functionality.
    f_hos = HomoSeqScorer(msa, f_5p, f_len, f_a, 100, False)
    r_hos = HomoSeqScorer(msa, r_5p, r_len, r_a, 100, True)

    f_hos.find_best()
    r_hos.find_best()

    assert len(f_hos.get_best_params()) == len(r_hos.get_best_params()) == 9

    hes = HeteroSeqScorer(msa, f_5p, r_5p, f_len, r_len, amp_lens, f_a, r_a, 100)

    hes.find_best()

    assert len(hes.get_best_params()) == 81


def test_complex_operation_0() -> None:
    """Tests that the program is successfully able to locate dimer
    structures. """
    msa = get_msa('25_seqs')
    amp_lens = list(range(60, 71))

    f_5p = list(range(380, 400))
    r_5p = list(range(440, 480))

    f_len = [15]
    r_len = [15]

    f_a = ['ACACTCTTTCCCTACACGACGCTCTTCCGATCT']
    r_a = ['GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT']

    # Basic HomoSeqScorer functionality.
    f_hos = HomoSeqScorer(msa, f_5p, f_len, f_a, 100, False)
    r_hos = HomoSeqScorer(msa, r_5p, r_len, r_a, 100, True)

    f_hos.find_best()
    r_hos.find_best()

    assert len(f_hos.get_best_params()) == 20
    assert len(r_hos.get_best_params()) == 40

    hes = HeteroSeqScorer(msa, f_5p, r_5p, f_len, r_len, amp_lens, f_a, r_a, 100)

    hes.find_best()

    assert len(hes.get_best_params()) == 100






