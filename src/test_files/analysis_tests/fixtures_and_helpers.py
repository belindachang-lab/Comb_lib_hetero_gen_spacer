from src.meta_tools.analysis_tools import Heterodimer, Homodimer
from src.meta_tools.evaluate import DimerManager
import random

# The attributes of a standard dimer manager

standard_fh = [
    Homodimer(0, True, 9, 12),
    Homodimer(0, True, -2, 22),
    Homodimer(0, True, 7, 4),
    Homodimer(1, True, -3, 2),
    Homodimer(1, True, 44, 12),
    Homodimer(1, True, -6, 11),
    Homodimer(2, True, 33, 15),
    Homodimer(2, True, -27, 4),
]

standard_rh = [
    Homodimer(2, True, -9, 12),
    Homodimer(0, True, 2, 22),
    Homodimer(0, True, -1, 4),
    Homodimer(3, True, 6, 2),
    Homodimer(1, True, -22, 12),
    Homodimer(1, True, 13, 11),
    Homodimer(2, True, -21, 15),
    Homodimer(3, True, 8, 4),
]

standard_he = [
    Heterodimer(2, 0, -30, 12),
    Heterodimer(0, 1, 4, 22),
    Heterodimer(0, 2, -29, 4),
    Heterodimer(3, 3, 30, 2),
    Heterodimer(1, 1, -65, 12),
    Heterodimer(1, 0, 55, 11),
    Heterodimer(2, 3, -2, 15),
    Heterodimer(3, 2, 7, 4),
]

STANDARD_DM = DimerManager(standard_fh, standard_rh, standard_he)


def get_rand_homo() -> Homodimer:
    """Returns a randomly generated homodimer. Homodimer will have allowable
    values."""
    ind = random.randrange(0, 3)
    pdir = random.randint(0, 1)
    off = random.randrange(0, 100)
    comp = random.randrange(-40, 40)
    return Homodimer(ind, bool(pdir), off, comp)


def get_rand_hetero() -> Heterodimer:
    """Returns a randomly generated homodimer. Homodimer will have allowable
    values."""
    f_ind = random.randrange(0, 3)
    r_ind = random.randrange(0, 3)
    off = random.randrange(0, 100)
    comp = random.randrange(-40, 40)
    return Heterodimer(f_ind, r_ind, off, comp)

rand_fh = []
rand_rh = []
rand_hd = []
num_rand = 1000

for q in range(num_rand):
    rand_fh.append(get_rand_homo())
    rand_rh.append(get_rand_homo())
    rand_hd.append(get_rand_hetero())
    rand_hd.append(get_rand_hetero())

RAND_DM = DimerManager(rand_fh, rand_rh, rand_hd)
