from seq_alignment_analyser import *
from src.test_files.fixtures_and_helpers import TEST_OUTPUT_PATH
from test_align import get_msa

to_run_name = '25_seqs'


msa = get_msa(to_run_name)
msa.show_plot()
msa.save_plot(TEST_OUTPUT_PATH / (to_run_name + '.png'))
