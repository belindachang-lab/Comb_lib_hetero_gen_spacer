from src.hetero_spacer_generator.spacer_generator import HeteroGen
from src.hetero_spacer_generator import eval_total_complementarity
from src.hetero_spacer_generator.spacer_generator.spacer_filters import \
    SortForSimultaneous, cross_compare, remove_high_consec_complementarity, \
    remove_high_dimer_complementarity
from src.test_files.fixtures_and_helpers import *
import pytest


class TestSortForSimultaneous:

    def test_remove_high_dimer_complementarity_basic(self) -> None:
        """Tests whether remove_high_dimer_complementarity produces a sample
        with the correct output size"""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        remove_high_dimer_complementarity(sfm.forward_spacer_seqs,
                                          sfm.incomplete_forward_primer,
                                          SMALL_SAMPLE_SIZE)
        remove_high_dimer_complementarity(sfm.reverse_spacer_seqs,
                                          sfm.incomplete_reverse_primer,
                                          SMALL_SAMPLE_SIZE)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_remove_high_dimer_complementarity_outprm_random(self) -> None:
        """Tests whether the remove_high_dimer_complementarity produces primers
        with less binding than a random sampling."""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12)
        rand_for_spacers = []
        for i in range(SMALL_SAMPLE_SIZE):
            rand_for_spacers.append(random.choice(sfm.forward_spacer_seqs))
        remove_high_dimer_complementarity(sfm.forward_spacer_seqs,
                                          sfm.incomplete_forward_primer,
                                          SMALL_SAMPLE_SIZE)
        rand_scores = 0
        sorted_scores = 0
        for i in range(SMALL_SAMPLE_SIZE):
            sorted_scores += eval_total_complementarity(
                sfm.incomplete_forward_primer,
                sfm.forward_spacer_seqs[i])
            rand_scores += eval_total_complementarity(
                sfm.incomplete_forward_primer,
                rand_for_spacers[i])
        assert rand_scores >= sorted_scores

    def test_remove_high_consec_complementarity_basic(self) -> None:
        """Tests whether remove_high_dimer_complementarity produces a sample with the
        correct output size"""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        remove_high_consec_complementarity(sfm.forward_spacer_seqs,
                                           sfm.incomplete_forward_primer,
                                           SMALL_SAMPLE_SIZE)
        remove_high_consec_complementarity(sfm.reverse_spacer_seqs,
                                           sfm.incomplete_reverse_primer,
                                           SMALL_SAMPLE_SIZE)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_filter_spacer_sets_basic(self) -> None:
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        sss = SortForSimultaneous(12, 12)
        sss._filter_spacer_sets(sfm.incomplete_forward_primer,
                                sfm.incomplete_reverse_primer,
                                sfm.forward_spacer_seqs,
                                sfm.reverse_spacer_seqs,
                                SMALL_SAMPLE_SIZE / MED_SAMPLE_SIZE * 100)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_cross_compare_basic(self):
        """Tests that cross compare prodcues the expected output for a basic
        input."""
        sfm = SeqFixtureManager()
        sfm.do_all()
        sss = SortForSimultaneous(12, 12)
        numsets = 3
        primer_sets = cross_compare(sfm.incomplete_forward_primer,
                                    sfm.incomplete_reverse_primer,
                                    sfm.forward_spacer_seqs[
                                    0:SMALL_SAMPLE_SIZE],
                                    sfm.reverse_spacer_seqs[
                                    0:SMALL_SAMPLE_SIZE],
                                    numsets)
        assert len(primer_sets) == numsets

    def test_full_primer_creation_process(self):
        """for_adapter_seq = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
        rev_adapter_seq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
        for_indexing_seq = 'ATCG'
        rev_indexing_seq = 'GCTA'
        for_binding_seq = 'TGTCATCTCCTTCTGCTACGG'
        rev_binding_seq = 'GCCGGCATGGTCATGAAG'"""
        for_adapter_seq = 'ACACTCWTTCCCTACACGACGNTCCGATCT'
        rev_adapter_seq = 'GTGACTGGAGTTNAGACGTGTGCTCTTCCGATCT'
        for_indexing_seq = 'AWCG'
        rev_indexing_seq = 'GBTA'
        for_binding_seq = 'TGTCATCVCCTTCTGCTACGG'
        rev_binding_seq = 'GCCGGCAVGGTCATGAAG'

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

        hg = HeteroGen(hetero_size, spacer_size)
        hg.set_rigour(-5)
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


if __name__ == '__main__':
    pytest.main(['test_files/simultaneous_tests.py', '-v'])
