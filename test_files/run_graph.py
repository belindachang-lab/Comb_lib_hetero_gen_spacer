from seq_alignment_analyser import *
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH, TEST_OUTPUT_PATH
from test_align import get_msa
import matplotlib.pyplot as plt

to_run_name = '25_seqs'


msa = get_msa(to_run_name)
msa.show_plot()
msa.save_plot(TEST_OUTPUT_PATH / (to_run_name + '.png'))
