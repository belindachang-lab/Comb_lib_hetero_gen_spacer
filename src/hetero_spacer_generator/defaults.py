from pathlib import Path

# === Misc ===

# Verbosity
from src.execution_managers.parameter_manager import get_pm

V = True

# Processes to spawn.
NUM_PROCS = 1

# Print runtime of methods.
TIMING = True

#  Store csv of primer output files.
SCORE_CSV = True

# The path to deposit the csv scores path.
CSV_PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\'
                'Scores CSV')

# === General Defaults ===
NUM_SPACERS = 4
# The default maximum length of any spacer created
MAX_SPACER_LENGTH = 12
ABSOLUTE_MAX_SPACER_LENGTH = 50
# The default size of the heterogeneity region to proceed the binding region of
# the primer
NUM_HETERO = 12
ABSOLUTE_MAX_NUM_HETERO = 50
RIGOUR = 0
# === Alignment Analyser Defaults ===
# Minimum complementarity at binding site in alignment.
MIN_COMP = 85
# Defualt window size when generating sliding window charts.
WINDOW_SIZE = 10
# The number of primers to keep during execution.
NUM_TO_KEEP = 10000

# === Random sampling defaults ===
# The number of forward and reverse primers to compare against each other.
NUM_PAIRINGS_TO_COMP = 35
# The number of primers in the initial sets of potential primers.
INITIAL_PRIMER_SET_SIZE = 1000

# === Primer evaluation defaults ===
# Relative importance of consecutive bases in a dimer.
CONSEC_WEIGHT = 3
# Relative inmportance of non-consecutive bases in a dimer.
TOTAL_WEIGHT = 1


# === Alignment analysis gen ===
# Increasing the below factors changes the importance exponentially.
# Relative importance of binding potential when evaluating a set of binding
# sequences.
DIMER_WEIGHT = 1
# Relative importance of conservation across selected region.
CONS_WEIGHT = 2
# If we have them, use user define attributes.
try:
    pm = get_pm()
    cw = pm.get('CONSERVATION_WEIGHT')
    dw = pm.get('DIMER_WEIGHT')
except ValueError:
    dw = DIMER_WEIGHT
    cw = CONS_WEIGHT

DIMER_WEIGHT = dw
CONS_WEIGHT = cw

# === SpacerAlignmentGen Specific defaults ===
# Begin default criteria weighting.
GET_SMALLEST_TOTAL_LEN_DEFAULT = 1
# We make this 4, the corresponding method is likely to return values a fourth
# of the size of get_smallest_total_len_list().
GET_SMALLEST_OF_ANY_SPACER_DEFAULT = 4
# Very important that the spacers aren't aligned.
GET_UNIQUE_SPACER_OFFSETS_DEFAULT = 4

# Degenerate primers are enabled.
DEGENERATE_PRIMERS = False

DEGEN_TO_POSSIBLE = {
    "A": "A",
    "C": "C",
    "T": "T",
    "G": "G",
    "U": "U",
    "W": "AT",
    "S": "CG",
    "M": "AC",
    "K": "GT",
    "R": "AG",
    "Y": "CT",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
    # Custom: To be ignored.
    "I": ""
}


BASES = ['A', 'T', 'C', 'G']

