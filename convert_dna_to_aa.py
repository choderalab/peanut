"""
Requires scikit-bio
https://github.com/biocore/scikit-bio
"""
from skbio.sequence import genetic_code

def dna_to_aa(sequence):
    """
    sequence (string) DNA sequence
    """
    orig_code = genetic_code(11)
    return orig_code.translate(sequence).sequence


def all_dna_point_mutants_to_aa(wt_sequence):
    """
    wt_sequence (string) DNA sequence
    """
    AA_sequences = set()
    orig_code = genetic_code(11)
    AA_sequences.add(orig_code.translate(wt_sequence).sequence)

    for k, char in enumerate(wt_sequence):
        for mutant in ['A','C','G','U']:
            if mutant == char:
                continue
            if k == 0:
                this_dna_string = ""+mutant+wt_sequence[k+1:]
            if k > 0 and k < len(wt_sequence):
                this_dna_string = ""+wt_sequence[0:k]+mutant+wt_sequence[k+1:]
            if k == len(wt_sequence):
                this_dna_string = ""+wt_sequence[0:k]+mutant
            this_sequence = orig_code.translate(this_dna_string).sequence
            if '*' in this_sequence:
                continue
            else:
                AA_sequences.add(this_sequence)
    return AA_sequences


