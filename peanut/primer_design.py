"""
Requires scikit-bio
https://github.com/biocore/scikit-bio
"""
from skbio.sequence import genetic_code
from skbio.sequence import DNASequence

def check_melting_temp(primer_sequence):
    N = len(primer_sequence)
    gc_percent = float(primer_sequence.count('g') + primer_sequence.count('c')) / N * 100.0
    mismatch_percent = 1.000 / N * 100.0
    melting_temp = 81.5 + 0.41*gc_percent - 675.0/N - mismatch_percent
    if melting_temp < 78.0:
        return False
    else:
        # should be actively dealing with this; for now just giving notification
        if primer_sequence[0] not in ['g','c'] or primer_sequence[-1] not in ['g','c']:
            print("Primer sequence should optimally terminate in one or more C or G bases.")
        if gc_percent < 40.0:
            print("GC out of range!")
            print(str(gc_percent)+"% GC")
        print("Melting temp: "+str(melting_temp)+"C")
        return True

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

    Desired mutation must require only a single nucleotide change
    """
    orig_code = genetic_code(11)
    sequence = sequence.upper()
    aa_sequence = orig_code.translate(sequence).sequence

    if not str(wt_res) == aa_sequence[res_num-first_res]:
        raise IOError("Desired residue not found -- check wildtype residue name and id, and first residue id")
    # start of codon of residue of interest is at (res_num - first_res)*3

    wt_codon = DNASequence(sequence[(res_num - first_res)*3:(res_num - first_res)*3+3])

    mut_codons = orig_code.synonyms[mut_res]
    mut_codon = None
    for codon in mut_codons:
        if wt_codon.distance(DNASequence(codon))*3 == 1:
            mut_codon = codon

    if not mut_codon:
        raise IOError("Cannot make desired mutant with a single base change")

    good_melting_temp = False
    surrounding_bases = 11
    while not good_melting_temp:
        start_ix = max(0,(res_num-first_res)*3-surrounding_bases)
        end_ix = min(len(sequence),(res_num+1-first_res)*3+surrounding_bases)
        if end_ix - start_ix > 45:
            print("Acceptable melting temp was not found")
            break
        forward_primer = sequence[start_ix:(res_num - first_res)*3]+mut_codon+sequence[(res_num+1 - first_res)*3:end_ix]
        forward_primer = forward_primer.lower()
        good_melting_temp = check_melting_temp(forward_primer)
        surrounding_bases+=1
    
    forward_sequence = DNASequence(forward_primer)
    reverse_sequence = forward_sequence.rc()
    
    reverse_primer = reverse_sequence.sequence

    return forward_primer, reverse_primer

def make_mutant(sequence,wt_res,res_num,mut_res,first_res=1):
    """
    sequence (string) DNA sequence
    wt_res (char) single letter amino acid code of wildtype residue to be mutated
    res_num (int) residue id number of residue to be mutated
    mut_res (char) single letter amino acid code of mutant residue
    first_res (int) residue id number of first residue in sequence (default = 1)

    DNA sequence needs to start with the first residue of the protein (no promoter, etc)
    take DNA sequence, convert to AA, define AA point mutant, find corresponding codon of wt and mut, output forward and reverse primers
    DNA sequence should be only the kinase domain
    """

    orig_code = genetic_code(11)
    sequence = sequence.upper()
    aa_sequence = orig_code.translate(sequence).sequence

    if not str(wt_res) == aa_sequence[res_num-first_res]:
        raise IOError("Desired residue not found -- check wildtype residue name and id, and first residue id")
    # start of codon of residue of interest is at (res_num - first_res)*3
    start_ix = max(0,(res_num-first_res-5)*3)
    end_ix = min(len(sequence),(res_num-first_res+5)*3)

    wt_codon = DNASequence(sequence[(res_num - first_res)*3:(res_num - first_res)*3+3])

    mut_codons = orig_code.synonyms[mut_res]
    mut_codons = [DNASequence(codon) for codon in mut_codons]
    distances = [wt_codon.distance(codon) for codon in mut_codons]

    # choose the codon that requires fewest changes
    mut_codon = mut_codons[distances.index(min(distances))].sequence

    forward_primer = sequence[start_ix:(res_num - first_res)*3]+mut_codon+sequence[(res_num+1 - first_res)*3:end_ix]
    forward_primer = forward_primer.lower()

    forward_sequence = DNASequence(forward_primer)
    reverse_sequence = forward_sequence.rc()

    reverse_primer = reverse_sequence.sequence

    return forward_primer, reverse_primer




