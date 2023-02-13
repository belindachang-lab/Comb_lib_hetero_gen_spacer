from src.execution_managers.parameter_manager import *
from src.test_files.fixtures_and_helpers import READ_MES_PATH

dummy_fasta_lines = [
    '@BEGIN_FASTA',
    '>1',
    'google ',
    'plex',
    '>2',
    'suff',
    'ering',
    '@END_FASTA'
]

def get_parameter_manager(filename: str) -> ParameterManager:
    """Retries the lines from the given readme file."""
    with open(READ_MES_PATH / (filename + '.txt'), 'r') as file:
        return ParameterManager(file.readlines())



def test_parse_fasta() -> None:
    lines = parse_fasta(dummy_fasta_lines)
    assert 'googleplex' in lines
    assert 'suffering' in lines

def test_list_int() -> None:
    s = '1, 2, 3, 4,5,6'

def test_fully_filled_in() -> None:
    """Tests the parameter manager's ability to parse a fully filled in read
    me file."""
    pm = get_parameter_manager('filled_in')
    assert pm.has_all

def test_empty() -> None:
    """Tests that the parameter manager can parse the default file."""
    pm = get_parameter_manager('empty')
    aa = pm.can_aa()
    sfbr = pm.can_sfbr()
    hsg = pm.can_hsg()
    assert not (aa or sfbr or hsg)

