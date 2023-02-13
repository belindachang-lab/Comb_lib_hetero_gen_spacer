import pytest
from seq_alignment_analyser.align import *
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH, get_msa
import os

all_alignment_paths = []
for test_file_path in os.listdir(ALIGNMENTS_PATH):
    all_alignment_paths.append(ALIGNMENTS_PATH / test_file_path)

def test_complementarity() -> None:
    """Tests correct functionality of complementarity function."""
    for base in 'ATCG':
        assert get_complementarity(base, base) == 1
        assert get_complementarity(base, 'N') == 1 / 4
        assert get_complementarity('N', base) == 1 / 4

    assert get_complementarity('Y', 'R') == 0

@pytest.mark.parametrize('align_path', all_alignment_paths)
def test_parse_adapter(align_path) -> None:
    """Tests that all alignments can be succesfully parsed into an MSA
    object."""
    msa = MSA(align_path)
    assert msa.get_num_seqs() > 1


def test_correct_seq_parsing() -> None:
    """Tests that the correct number of sequences are parsed."""
    msa = get_msa('25_seqs')
    assert msa.get_num_seqs() == 25

@pytest.fixture
def cons_score_msa() -> MSA:
    return get_msa('conservation')

def test_full_conservation(cons_score_msa) -> None:
    """Tests the proper calculation of full conservation scores."""
    assert cons_score_msa.get_conservation(0) == 100


def test_partial_conservation(cons_score_msa) -> None:
    """Tests the proper calculation of partial conservation scores."""
    assert cons_score_msa.get_conservation(1) == 25
    assert cons_score_msa.get_conservation(2) == 25

def test_prcnt_spacer(cons_score_msa) -> None:
    """Tests proper calculation of percent spacer scores."""
    assert cons_score_msa.get_percent_spacer(2) == 75
    assert cons_score_msa.get_percent_spacer(0) == 0

def test_percent_missed(cons_score_msa) -> None:
    """Tests proper calculation of the percent of bases missed."""
    assert cons_score_msa.get_percent_missed(0) == 0
    assert cons_score_msa.get_percent_missed(2) == 25
    assert cons_score_msa.get_percent_missed(4) == 50
    assert cons_score_msa.get_percent_missed(5) == 25

def test_warning_str(cons_score_msa) -> None:
    s = cons_score_msa.scan_region(0, 10)
    print(s)
    assert s
    assert False
