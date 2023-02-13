from src.hetero_spacer_generator.sequence_tools import is_degen
from src.execution_managers import demo_tools
from time import time



hg = demo_tools.HeteroGen()
input("Press enter to start.")

demo_tools.ACCEPT_DEGEN = True


incomplete_forward_primer, forward_spacer = demo_tools.get_primer_and_spacers(hg,
                                                                   'forward',
                                                                              is_one_step=True)
incomplete_reverse_primer, reverse_spacer = demo_tools.get_primer_and_spacers(hg,
                                                                   'reverse',
                                                                              is_one_step=True)

while True:

    degen = False
    for primer in (incomplete_forward_primer, incomplete_reverse_primer):
        for seq in primer:
            if is_degen(seq):
                degen = True
                break

    hg.set_pairwise()

    rigour = int(demo_tools.while_not_valid(
        "Enter the rigour with which the program "
        "should search for possible combinations of spacers. "
        "Note: increasing the rigour leads to "
        "an exponential increase in runtime. Enter a value between -10 and 10:",
        "Please enter a number between -10 and 10:",
        demo_tools.RANGE, start=-10, end=10))

    hg.set_rigour(rigour)

    number_to_return = int(
        demo_tools.while_not_valid("Enter the number of sets of primer sets you'd like to "
                        "see.", "Please enter a number between 1 and 10:",
                                   demo_tools.RANGE, start=1, end=10))

    print("Generating primers. This may take some time...")

    t0 = time()

    psets = hg.get_hetero_seqs(incomplete_forward_primer, incomplete_reverse_primer,
                       forward_spacer, reverse_spacer, number_to_return)

    runtime = time() - t0

    input("Completed search in {time:.2f} seconds. Press enter t"
          "o see results.".format(time=runtime))

    for i in range(len(psets)):
        print("Set #", i + 1, ':')
        print(str(psets[i]))

    input("Press enter to generate a new set of primers.")






