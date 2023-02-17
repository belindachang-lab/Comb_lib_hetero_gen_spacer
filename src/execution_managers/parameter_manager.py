from pathlib import Path
from typing import List, Any, Dict, Tuple
import os

class EmptyContents(Exception):
    """Raised when the contents of a variable are empty."""
    pass

def parse_int(s: str) -> int:
    """Parses a string from the config_read_me, converting it into the required
    type."""
    return int(s)


def parse_str(s: str) -> str:
    """Parses a string from the config_read_me, converting it into the required
    type."""
    return s.strip()


def contains(d1: dict, d2: dict) -> bool:
    """Returns whether all keys in d2 are in d1."""
    for key in d2.keys():
        if key not in d1.keys():
            return False
    return True

def get_not_in(d1: dict, d2: dict) -> List[str]:
    """Returns whether all keys in d2 are in d1."""
    not_in = []
    for key in d2.keys():
        if key not in d1.keys():
            not_in.append(key)
    return not_in


def parse_bool(s: str) -> bool:
    """Parses a string from the config_read_me, converting it into the required
    type."""
    if 'T' in s.upper():
        return True
    elif 'F' in s.upper():
        return False
    else:
        raise ValueError('')


def parse_fasta(lines: List[str]) -> List[str]:
    """Extracts a fasta from the given lines. Assumes that lines begins with
    @BEGIN_FASTA and ends with @END_FASTA."""
    # Error checking...
    if not len(lines) > 2:
        raise EmptyContents('')
    is_begin = "@BEGIN_FASTA" in lines[0]
    is_end = "@END_FASTA" in lines[-1]
    proper_start = '>' in lines[1]

    if not (is_begin and is_end and proper_start):
        raise ValueError('')

    # Remove begin and end lines.
    lines.pop(0)
    lines.pop(-1)

    # Extract strings.
    strz = []
    cur_s = ''
    for line in lines:
        if '>' in line:
            strz.append(cur_s)
            cur_s = ''
            continue
        cur_s += line.strip()
    strz.append(cur_s)
    strz.pop(0)

    if len(strz) == 0:
        raise EmptyContents('')

    return strz


def parse_list_int(s: str) -> List[int]:
    """Parses a string from the config_read_me, converting it into the required
    type."""
    ints = s.split(',')
    ints = [n.strip() for n in ints]
    ints = [int(n) for n in ints]
    set_ints = set(ints)
    if len(set_ints) != len(ints):
        raise ValueError('Duplicate integers.')
    return ints


def filter_lines(lines: List[str]) -> None:
    """Removes unwanted lines from <lines>."""
    for i in range(len(lines) - 1, -1, -1):
        line = lines[i].strip()
        if line == '':
            lines.pop(i)
        elif line[0] == '#':
            lines.pop(i)
        elif line[0] == '*':
            lines[i] = line[1:]
        else:
            lines[i] = line


def fetch_fasta_lines(lines: List[str], start: int) -> Tuple[List[str], int]:
    """Extracts fasta from lines, where lines[start] contains @BEGIN_FASTA.
    Returns extracted lines as well as finish lines containing @END_FASTA."""
    i = start

    while not '@BEGIN_FASTA' in lines[i]:
        i += 1
    fasta = []
    fasta.append(lines[i])
    while not '@END_FASTA' in lines[i]:
        i += 1
        fasta.append(lines[i])

    return fasta, i


class ParameterManager:
    """A class designed to parse parameters in a textfile and return them when
    needed throughout the program.

    === Private Attributes ===

    _values_dict:
        Stores all values parsed from the config readme.
    _pot_values_dict:
        Stores all possible values and their types.

    _req_aa:
        All parameters required to run analyse alignment.

    _req_sfbr:
        All parameters required to run scan for binding regions.

    _req_hsg:
        All parameters required to run heterogeneity sequence generator.
    """

    _values_dict: Dict[str, Any]
    _pot_values_dict: Dict[str, str]
    _req_aa: Dict[str, str]
    _req_sfbr: Dict[str, str]
    _req_hsg: Dict[str, str]

    _lines: List[str]
    _i: int

    def __init__(self, lines: List[str]) -> None:
        """Parses the provided lines into a dictionary of values."""
        self._init_reqs()
        self._parse_readme_lines(lines)

    def has(self, var_name) -> bool:
        """Returns whether a value for <var_name> has been received."""
        return var_name in self._values_dict.keys()

    def get(self, var_name: str):
        """Gets the value of the requested variable."""
        try:
            return self._values_dict[var_name]
        except KeyError:
            raise ValueError('Value for ' + str(var_name) +
                             ' has not been received.')

    def _parse_readme_lines(self, lines: List[str]) -> None:
        """Converts the values stored in lines into the values dict."""
        filter_lines(lines)
        self._lines = lines
        self._values_dict = {}
        self._i = 0
        while self._i < len(lines):
            line = lines[self._i]
            var_name = line.split('=')[0]
            var_name = var_name.strip()
            self._parse_line(line, var_name)
            self._i += 1

        self._misc()

    def _misc(self) -> None:
        """Perform minor adjustments to the received data."""

        if self.get('OUTPUT_FASTA') == 'DIR':
            self._values_dict['OUTPUT_FASTA'] = \
                Path(os.path.dirname(__file__)).parent.parent / \
                'output_primers.fasta'

        if self.get('NUM_CORES') == 0:
            self._values_dict['NUM_CORES'] = os.cpu_count()

    def has_all(self) -> None:
        """Checks that all required values are present."""
        var_names = list(self._pot_values_dict.keys())
        for key in self._values_dict.keys():
            var_names.remove(key)

        if len(var_names) > 0:
            raise ValueError('Not all values received. Missing: ' +
                             str(var_names))


    def _parse_line(self, line: str, var_name: str) -> None:
        """Inserts the information contained in the given line into
        values_dict."""

        contents = line.split('=')[-1]

        # No value specified.
        if contents == var_name or '#' in line:
            return

        try:
            var_type = self._pot_values_dict[var_name]
        except KeyError:
            raise ValueError('Invalid input, offending line: \'' + line + '\'')

        match var_type:
            case 'i':
                parse_fxn = parse_int
            case 's':
                parse_fxn = parse_str
            case 'b':
                parse_fxn = parse_bool
            case 'ls':
                contents, self._i = fetch_fasta_lines(self._lines, self._i)
                parse_fxn = parse_fasta
            case 'li':
                parse_fxn = parse_list_int
            case _:
                raise ValueError('I have no idea how you did this.')

        try:
            self._values_dict[var_name] = parse_fxn(contents)
        except EmptyContents:
            pass
        except ValueError:
            raise ValueError('Invalid input, offending line: \'' + line + '\'')

    def can_aa(self) -> bool:
        """Returns whether the variables neccesary to run analyse alignment
         have been received."""
        return contains(self._values_dict, self._req_aa)

    def can_sfbr(self) -> bool:
        """Returns whether the variables neccesary to run scan for binding
        regions have been received."""
        return contains(self._values_dict, self._req_sfbr)

    def can_hsg(self) -> bool:
        """Returns whether the variables neccesary to run heterogeneity
         spacer gen have been received."""
        return contains(self._values_dict, self._req_hsg)

    def missing_aa(self) -> List[str]:
        """Returns whether the variables neccesary to run analyse alignment
         have been received."""
        return get_not_in(self._values_dict, self._req_aa)

    def missing_sfbr(self) -> List[str]:
        """Returns whether the variables neccesary to run scan for binding
        regions have been received."""
        return get_not_in(self._values_dict, self._req_sfbr)

    def missing_hsg(self) -> List[str]:
        """Returns whether the variables neccesary to run heterogeneity
         spacer gen have been received."""
        return get_not_in(self._values_dict, self._req_hsg)

    def _init_reqs(self) -> None:
        """Initialises each of the requirement dictionaries. Each dict maps
        name of variable to type, where
        i: int
        s: str
        b: bool
        ls: List[str]
        li: List[int]
        Prefacing with 'o' denotes the value is optional.
        """
        self._req_aa = {
            'MSA_FILEPATH': 's',
            'GRAPH_OUT_FILEPATH': 's',
            'STORE_IMAGE': 'b',
            'DISPLAY_GRAPH': 'b',
            'WINDOW_SIZE': 'i'
        }

        self._req_sfbr = {
            'MSA_FILEPATH': 's',
            'FORWARD_ADAPTERS': 'ls',
            'REVERSE_ADAPTERS': 'ls',
            'FORWARD_BINDING_START': 'i',
            'FORWARD_BINDING_END': 'i',
            'REVERSE_BINDING_START': 'i',
            'REVERSE_BINDING_END': 'i',
            'FORWARD_BINDING_LENGTHS': 'li',
            'REVERSE_BINDING_LENGTHS': 'li',
            'AMPLICON_LENGTH_MIN': 'i',
            'AMPLICON_LENGTH_MAX': 'i',
            'CONSERVATION_WEIGHT': 'i',
            'DIMER_WEIGHT': 'i',
        }

        self._req_hsg = {
            'FORWARD_ADAPTER': 's',
            'REVERSE_ADAPTER': 's',
            'FORWARD_BINDING_SEQ': 's',
            'REVERSE_BINDING_SEQ': 's',
            'HETEROGENEITY_REGION_LENGTH': 'i',
            'OUTPUT_FASTA': 's',
            'SHOW_SPACER_MENU': 'b',
            'NUM_CORES': 'i',
            'RIGOUR': 'i',
            'NUM_SETS_TO_GEN' : 'i'
        }

        self._pot_values_dict = {}
        for dic in (self._req_aa, self._req_sfbr, self._req_hsg):
            for key in dic.keys():
                self._pot_values_dict[key] = dic[key]

def get_pm() -> ParameterManager:
    """Generates a parameter manager by parsing the readme in the directory."""
    dir_path = Path(os.path.dirname(__file__)).parent.parent
    read_me_path = dir_path / 'config_read_me.txt'
    with open(read_me_path, 'r') as file:
        return ParameterManager(file.readlines())

