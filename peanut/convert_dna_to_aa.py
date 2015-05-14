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
    rv_sequence = sequence[::-1]

    translated = []
    stops = []

    for i in range(3):
        translated.append(orig_code.translate(sequence).sequence)
        stops.append(translated[-1].count('*'))
        translated.append(orig_code.translate(rv_sequence).sequence)
        stops.append(translated[-1].count('*'))
        sequence = sequence[1:]
        rv_sequence = rv_sequence[1:]

    return translated[stops.index(min(stops))]


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

