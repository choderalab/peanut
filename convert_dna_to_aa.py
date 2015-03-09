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


def two_dna_point_mutants_to_aa(wt_sequence):
    """
    wt_sequence (string) DNA sequence
    """
    AA_sequences = set()
    orig_code = genetic_code(11)
    AA_sequences.add(orig_code.translate(wt_sequence).sequence)

    for k1, char1 in enumerate(wt_sequence[:-1]):
        for k2, char2 in enumerate(wt_sequence[k1+1:]):
            for mutant1 in ['A','C','G','U']:
                if mutant1 == char1:
                    continue
                for mutant2 in ['A','C','G','U']:
                    if k1 == 0:
                        if k2 < len(wt_sequence):
                            this_dna_string = ""+mutant1+wt_sequence[k1+1:k2]+mutant2+wt_sequence[k2+1:]
                        else:
                            this_dna_string = ""+mutant1+wt_sequence[k1+1:k2]+mutant2
                    if k1 > 0 and k1 < len(wt_sequence)-1:
                        if k2 < len(wt_sequence):
                            this_dna_string = ""+wt_sequence[0:k1]+mutant1+wt_sequence[k1+1:k2]+mutant2+wt_sequence[k2+1:]
                        else:
                            this_dna_string = ""+wt_sequence[0:k1]+mutant1+wt_sequence[k1+1:k2]+mutant2
                    if k1 == len(wt_sequence)-1:
                        this_dna_string = ""+wt_sequence[0:k1]+mutant1+mutant2
                    this_sequence = orig_code.translate(this_dna_string).sequence
                    if '*' in this_sequence:
                        continue
                    else:
                        AA_sequences.add(this_sequence)
    return AA_sequences

def calculate_extinction_coeff(aa_sequence):
    """
    Extinction Coefficient in cm^-1 M^-1
    """
    count_Tyr = aa_sequence.count('Y')
    count_Trp = aa_sequence.count('W')
    count_Cys = aa_sequence.count('C')
    cytosolic = 1490*count_Tyr + 5500*count_Trp
    secreted = 1490*count_Tyr + 5500*count_Trp + 125*count_Cys/2
    return secreted, cytosolic

def calculate_molecular_weight(aa_sequence):
    """
    Molecular Weight in g/mol
    Note that this doesn't agree with that website
    """
    # wtf do you do with B and Z?
    aa_to_mw = {'A':89.0935, 'C':121.1590, 'D':133.1032, 'E':147.1299, 'F':165.1900, 'G':75.0669, 'H':155.1552, 'I':131.1736, 'K':146.1882, 'L':131.1736, 'M':149.2124, 'N':132.1184, 'P':115.1310, 'Q':146.1451, 'R':174.2017, 'S':105.0930, 'T':119.1197, 'V':117.1469, 'W':204.2262, 'Y':181.1894}
    molecular_weight = 0.0
    for k, char in enumerate(aa_sequence):
        molecular_weight += aa_to_mw[char]
        if k != 0: # water mlc given off in peptide bond
            molecular_weight -= 18.01528
    return molecular_weight

def calculate_absorbance(aa_sequence):
    secreted, cytosolic = calculate_extinction_coeff(aa_sequence)
    molecular_weight = calculate_molecular_weight(aa_sequence)
    # which extinction coeff do you use?
    # or both?
    absorbance = secreted / molecular_weight
    return absorbance

