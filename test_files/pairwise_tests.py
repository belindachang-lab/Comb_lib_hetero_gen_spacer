import random

import pytest
import fixtures_and_helpers as fah
from src.hetero_spacer_generator import MBPrimer, \
    MBPrimerBuilder, PairwisePrimerSet
from src.hetero_spacer_generator.spacer_generator.criteria import EvalMBPrimer
from src.hetero_spacer_generator.sequence_tools import SeqAnalyzer
from src.hetero_spacer_generator.spacer_generator import \
    HeteroGen, SpacerAlignmentGen
from src.hetero_spacer_generator.spacer_generator.random_spacer_generator import \
    gen_hetero_set
from src.hetero_spacer_generator.spacer_generator.spacer_filters import \
    SortForPairwise
from src.hetero_spacer_generator.spacer_generator.random_spacer_generator import \
    RandomSpacerGen

from src.meta_tools.template_sequences import *

class TestPairwiseSorter:
    sfm = fah.SeqFixtureManager()
    sfm.num_to_generate = fah.MED_SAMPLE_SIZE
    sfm.do_all()


DEF_RIGOUR = -5

@pytest.mark.parametrize('gen_set',
                         [
                             (STANDARD_SET),
                             (RANDOM_RHO_SET),
                             (SMALL_BINDING_SET),
                             (OPTIMAL_SET)
                         ], scope='class')
def test_increased_performace_rand_seq_sel(gen_set) -> None:
    """Tests that primers produced by the full filtering process
    consistently from less dimers than randomly generated primers of the
    same variety. Does not vary the heterogeneity combo used - will allways
    use the best available one."""
    num_to_test = 1
    f_bind, r_bind, f_adapt, r_adapt = gen_set.unpack()

    f_incomp = MBPrimerBuilder(binding_seq=f_bind[0],
                               adapter_seq=f_adapt[0])
    r_incomp = MBPrimerBuilder(binding_seq=r_bind[0],
                               adapter_seq=r_adapt[0])

    # Default length, something else otherwise.
    heterogeneity_region_size = min(12, min(len(f_bind[0]) - 1,
                                            len(r_bind[0])) - 1)
    rsg = RandomSpacerGen(heterogeneity_region_size,
                          heterogeneity_region_size, rigour=DEF_RIGOUR)
    sag = SpacerAlignmentGen(heterogeneity_region_size,
                             heterogeneity_region_size)
    f_spacer = sag.get_all_spacer_combos(f_bind[0])[0]
    r_spacer = sag.get_all_spacer_combos(r_bind[0])[0]
    for i in range(num_to_test):
        rand = rsg.get_random_hetero_seqs(f_incomp, r_incomp, f_spacer,
                                          r_spacer, 1)[0]
        filtered = rsg.get_hetero_seqs(f_incomp, r_incomp, f_spacer,
                                       r_spacer, 1)[0]

        rand_score = rand.get_score()
        filtered_score = filtered.get_score()
        # The score is proportional to the most stable heterodimer.
        # Higher => Worse
        assert rand_score > filtered_score

class TestSeqAnalyser:

    sa = SeqAnalyzer(degen=False)


    def test_com_seqs_basic(self) -> None:
        """Tests that sequence comparison works on basic sequences."""
        s1 = Seq('AAAAA')
        s2 = Seq('AAAAA')
        assert self.sa.comp_seqs_any_overlap(s1, s2, sa.get_consec_complementarity) == 0

        s1 = Seq('AAAAA')
        s2 = Seq('TTTTT')
        assert self.sa.comp_seqs_any_overlap(s1, s2, sa.get_consec_complementarity) == 5

        s1 = Seq('AATT')
        s2 = Seq('AATT')
        assert self.sa.comp_seqs_any_overlap(s1, s2, sa.get_consec_complementarity) == 4


    def test_com_seqs_any_overlap_prod_seqs_same_len(self) -> None:

        def equal_length(seq1: str, seq2: str) -> int:
            assert len(seq1) == len(seq2)
            return 1

        seqa = SeqAnalyzer()

        for i in range(fah.SMALL_SAMPLE_SIZE):
            seq1 = fah.gen_random_seq(random.randrange(0, fah.SMALL_SEQ_LEN))
            seq2 = fah.gen_random_seq(random.randrange(0, fah.SMALL_SEQ_LEN))
            seqa.comp_seqs_any_overlap(seq1, seq2[::-1], equal_length)

    def test_consec_complementarity(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(False)

        seq1 = Seq('ACAAGTGCTAGTAGCTGATCG')  # 5' - 3'
        seq2 = Seq('TTGATGAGATCAATGCTGTAG')  # 3' - 5'
        # Matching         ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

    def test_consec_complementarity_ends(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(False)

        seq1 = Seq('CTAGTACAAGTGAGCTGATCG')  # 5' - 3'
        # Matching  ^^^^^
        seq2 = Seq('TTGATGAATGCTGTAGGATCA')  # 3' - 5'
        # Matching                  ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

        seq1 = Seq('ACAAGTGAGCTGATCGCTAGT')  # 5' - 3'
        # Matching                  ^^^^^
        seq2 = Seq('GATCATTGATGAATGCTGTAG')  # 3' - 5'
        # Matching  ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

    def test_consec_complementarity_degen(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(True)

        seq1 = Seq('AKAAGTGCTNGTAGCTGATNG')  # 5' - 3'
        seq2 = Seq('TTGATGAGANCAATGCTGTNG')  # 3' - 5'
        # Matching         ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

    def test_consec_complementarity_ends_degen(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(True)

        seq1 = Seq('CTAGTACAAGTGAGCTGATCG')  # 5' - 3'
        # Matching  ^^^^^
        seq2 = Seq('TTGATGAATGCTGTAGGATCA')  # 3' - 5'
        # Matching                  ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

        seq1 = Seq('ACAAGTGAGCTGATCGCTVGT')  # 5' - 3'
        # Matching                  ^^^^^
        seq2 = Seq('GATCRTTGATGAATGCTGTAG')  # 3' - 5'
        # Matching  ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity) == 5

    def test_consec_complementarity_ends_with_skip(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(False)

        seq2 = Seq('ATATTGAATGCTGTAGGATCA')  # 3' - 5'
        # Matching  ^^^^4           ^^^^^5
        seq1 = Seq('CTAGTACAAGTGAGCTGTATA')  # 5' - 3'
        # Matching  ^^^^^5           ^^^^4

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity, 5,
                                          0) == 4

        assert seqa.comp_seqs_any_overlap(seq1, seq2[::-1],
                                          seqa.get_consec_complementarity, 0,
                                          5) == 5


# Primer regions with predetermined binding for testing
ads = Seq("TAGAGAGAGAATG")
ids = Seq("ATCG")
bin = Seq("ACTATATATATCG")
#            ^^^^^^^^^ 8
b1 = MBPrimerBuilder(ads, ids, binding_seq=bin)

# Primer regions with predetermined binding for testing
ads = Seq("ATCTCTCTCAATCGA")
#           ^^^^^^^^ 8
ids = Seq("ATCG")
bin = Seq("AAGAGAGAGAAGATTT")
#            ^^^^^^^^ 8
b2 = MBPrimerBuilder(ads, ids, binding_seq=bin)

ads = "TATATATATA"
ids = ""
hgs = "GGGG"
bin = "GCGCGCGCGC"

b3 = MBPrimer(ads, ids, hgs, bin)

ads = "TATAAAAAAA"
ids = ''
hgs = "GGGG"
bin = "AAAAAAGCGC"

b4 = MBPrimer(ads, ids, hgs, bin)


class TestEvalMBPrimer:
    """Test suite for EvalMBPrimer class"""
    emp = EvalMBPrimer(12, 12)

    def test_eval_inherent_homodimer_consec(self) -> None:
        assert self.emp.eval_inherent_homodimer_consec(b1) == 8
        assert self.emp.eval_inherent_homodimer_consec(b2) == 8

    def test_eval_inherent_heterodimer_consec(self) -> None:
        assert self.emp.eval_inherent_heterodimer_consec(b1, b2) == 8

    def test_eval_hetero_hetero_spacer_binding(self) -> None:
        assert self.emp.eval_homo_hetero_spacer_binding_consec(b3) == 8
        assert self.emp.eval_homo_hetero_spacer_binding_consec(b4) == 1


for_adapter_seq = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
rev_adapter_seq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
for_indexing_seq = 'ATCG'
rev_indexing_seq = 'GCTA'
for_binding_seq = 'TGTCATCTCCTTCTGCTACGG'
rev_binding_seq = 'GCCGGCATGGTCATGAAG'

hetero_size = 12
spacer_size = 12

incomp_forward_primer = MBPrimerBuilder()
incomp_forward_primer.set_binding_seq(for_binding_seq)
incomp_forward_primer.set_index_seq(for_indexing_seq)
incomp_forward_primer.set_adapter_seq(for_adapter_seq)
incomp_reverse_primer = MBPrimerBuilder()
incomp_reverse_primer.set_binding_seq(rev_binding_seq)
incomp_reverse_primer.set_index_seq(rev_indexing_seq)
incomp_reverse_primer.set_adapter_seq(rev_adapter_seq)


class TestSortForPairWise:

    def test_single_primer_sorting(self) -> None:
        hg = HeteroGen(hetero_size, spacer_size, rigour=DEF_RIGOUR)
        for_spacers = hg.get_all_spacer_combos(incomp_forward_primer
                                               .get_binding_seq())
        rev_spacers = hg.get_all_spacer_combos(incomp_reverse_primer
                                               .get_binding_seq())
        forward_spacer_seqs = gen_hetero_set(incomp_forward_primer,
                                             for_spacers[0],
                                             fah.LARGE_SAMPLE_SIZE)
        reverse_spacer_seqs = gen_hetero_set(incomp_reverse_primer,
                                             rev_spacers[0],
                                             fah.LARGE_SAMPLE_SIZE)
        ss = SortForPairwise(hetero_size, spacer_size, 10, False)
        ss._build_full(12, 12, incomp_forward_primer, incomp_reverse_primer,
                       forward_spacer_seqs, reverse_spacer_seqs)
        ss._evaluate_scores_single()
        ss._sort_and_trim()
        assert len(ss._for_seqs) == len(ss._rev_seqs) \
               == ss._num_pairings_to_compare

    def test_primer_set_sorting(self) -> None:
        pass

    def test_full_primer_creation_process(self):
        for_adapter_seq = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
        rev_adapter_seq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
        for_indexing_seq = 'ATCG'
        rev_indexing_seq = 'GCTA'
        for_binding_seq = 'TGTCATCTCCTTCTGCTACGG'
        rev_binding_seq = 'GCCGGCATGGTCATGAAG'

        hetero_size = 12
        spacer_size = 12

        incomp_forward_primer = MBPrimerBuilder()
        incomp_forward_primer.set_binding_seq(for_binding_seq)
        incomp_forward_primer.set_index_seq(for_indexing_seq)
        incomp_forward_primer.set_adapter_seq(for_adapter_seq)
        incomp_reverse_primer = MBPrimerBuilder()
        incomp_reverse_primer.set_binding_seq(rev_binding_seq)
        incomp_reverse_primer.set_index_seq(rev_indexing_seq)
        incomp_reverse_primer.set_adapter_seq(rev_adapter_seq)
        hg = HeteroGen(hetero_size, spacer_size, rigour=DEF_RIGOUR)
        for_spacers = hg.get_all_spacer_combos(incomp_forward_primer
                                               .get_binding_seq())
        rev_spacers = hg.get_all_spacer_combos(incomp_reverse_primer
                                               .get_binding_seq())

        primers = hg.get_hetero_seqs(incomp_forward_primer,
                                     incomp_reverse_primer,
                                     for_spacers[1],
                                     rev_spacers[1],
                                     5)
        assert len(primers) == 5
        for primer_set in primers:
            assert isinstance(primer_set, PairwisePrimerSet)









