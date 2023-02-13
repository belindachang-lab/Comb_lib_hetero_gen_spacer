from multiprocessing import freeze_support

from src.execution_managers.demo_tools import get_spacer
from src.execution_managers.parameter_manager import get_pm
from src.hetero_spacer_generator.primer_tools import MBPrimerBuilder
from src.hetero_spacer_generator.spacer_generator.hetero_spacer_generator import \
    HeteroGen
from src.meta_tools.gen_n_primers import primer_sets_to_str




def main():

    pm = get_pm()

    # Ensure that necessary parameters have been defined.
    if not pm.can_hsg():
        raise ValueError('Ensure all required parameters in the config file are '
                         'properly defined. Currently missing values include:' +
                         str(pm.missing_hsg())[1:-1])

    FORWARD_ADAPTER = pm.get('FORWARD_ADAPTER')
    REVERSE_ADAPTER = pm.get('REVERSE_ADAPTER')
    FORWARD_BINDING_SEQ = pm.get("FORWARD_BINDING_SEQ")
    REVERSE_BINDING_SEQ = pm.get("REVERSE_BINDING_SEQ")
    HETEROGENEITY_REGION_LENGTH = pm.get("HETEROGENEITY_REGION_LENGTH")
    OUTPUT_FASTA = pm.get("OUTPUT_FASTA")
    SHOW_SPACER_MENU = pm.get("SHOW_SPACER_MENU")
    NUM_CORES = pm.get("NUM_CORES")
    RIGOUR = pm.get("RIGOUR")
    NUM_SETS_TO_GEN = pm.get('NUM_SETS_TO_GEN')
    DIMER_WEIGHT = 1
    CONS_WEIGHT = 2

    # Get our spacer generator.
    hg = HeteroGen(rigour=RIGOUR,
                   num_hetero=HETEROGENEITY_REGION_LENGTH,
                   max_spacer_length=HETEROGENEITY_REGION_LENGTH)
    rg = hg.get_primer_gen()

    # Form incomplete primers. To be completed with heterogeneity spacer.
    incomp_f_primer = MBPrimerBuilder(adapter_seq=FORWARD_ADAPTER,
                                      binding_seq=FORWARD_BINDING_SEQ)
    incomp_r_primer = MBPrimerBuilder(adapter_seq=REVERSE_ADAPTER,
                                      binding_seq=REVERSE_BINDING_SEQ)

    # Obtain the sizing of the spacers.
    f_spacer = get_spacer(hg, incomp_f_primer, not SHOW_SPACER_MENU)
    r_spacer = get_spacer(hg, incomp_r_primer, not SHOW_SPACER_MENU)

    # Generate Primer Sets
    psets = rg.get_hetero_seqs(incomp_f_primer, incomp_r_primer,
                               f_spacer, r_spacer, NUM_SETS_TO_GEN, NUM_CORES)

    # Output to fasta.
    fasta_text = primer_sets_to_str(psets, True, 0, '\n')
    with open(OUTPUT_FASTA, "w") as my_file:
        my_file.write(fasta_text)

    print(NUM_SETS_TO_GEN, 'sets of primers output to \'' +
          str(OUTPUT_FASTA) + '\'')

if __name__ == '__main__':
    freeze_support()
    main()
