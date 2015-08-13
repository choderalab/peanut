"""
Requires scikit-bio
https://github.com/biocore/scikit-bio
"""
from skbio.sequence import genetic_code

def make_single_mutant(sequence,wt_res,res_num,mut_res,first_res=1):
    """
    sequence (string) DNA sequence
    wt_res (char) single letter amino acid code of wildtype residue to be mutated
    res_num (int) residue id number of residue to be mutated
    mut_res (char) single letter amino acid code of mutant residue
    first_res (int) residue id number of first residue in sequence (default = 1)

    DNA sequence needs to start with the first residue of the protein (no promoter, etc)
    take DNA sequence, convert to AA, define AA point mutant, find corresponding codon of wt and mut, output forward and reverse primers
    DNA sequence should be only the kinase domain
    need a way to define the AA numbering scheme / identify the correct residue
    """
    orig_code = genetic_code(11)
    aa_sequence = orig_code.translate(sequence).sequence

    if not str(wt_res) == aa_sequence[res_num-first_res]:
        raise IOError("Desired residue not found -- check wildtype residue name and id, and first residue id")
    # start of codon of residue of interest is at (res_num - first_res)*3
    start_ix = max(0,(res_num-first_res-5)*3)
    end_ix = min(len(sequence),(res_num-first_res+5)*3)
    wt_forward_primer = sequence[start_ix:end_ix]

    wt_codon = sequence[(res_num - first_res)*3:(res_num - first_res)*3+3]

    # there MUST be an actual function in scikit-bio that already does this
    mut_codon = None
    for k in range(3):
        for mutant in ['a','c','g','t']:
            if k == 0:
                new_codon = mutant+wt_codon[1:]
            if k == 1:
                new_codon = wt_codon[0]+mutant+wt_codon[2]
            if k == 2:
                new_codon = wt_codon[:-1]+mutant
            if orig_code.translate(new_codon).sequence == mut_res:
                mut_codon = new_codon
                break

    if not mut_codon:
        raise IOError("Cannot make desired mutant with a single base change")
    
    forward_primer = sequence[start_ix:(res_num - first_res)*3]+mut_codon+sequence[(res_num+1 - first_res)*3:end_ix]
    # HOW TO MAKE REVERSE PRIMER
    
    reverse_primer = "PLACEHOLDER" # PLACEHOLDER

    return forward_primer, reverse_primer
