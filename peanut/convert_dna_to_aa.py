"""
Requires scikit-bio
https://github.com/biocore/scikit-bio
"""
from skbio.sequence import genetic_code
from skbio.sequence import DNASequence

class SequenceConverter(object):
    def __init__(self, wt_dna_sequence):
        self.wt_dna_sequence = wt_dna_sequence

    def dna_to_aa(self, try_frames=False):
        aa_sequence = dna_to_aa(self.wt_dna_sequence, try_frames=try_frames)
        return aa_sequence
    
    def all_dna_point_mutants_to_aa(self):
        aa_sequences = all_dna_point_mutants_to_aa(self.wt_dna_sequence)
        return aa_sequences

    def two_dna_point_mutants_to_aa(self):
        aa_sequences = two_dna_point_mutants_to_aa(self.wt_dna_sequence)
        return aa_sequences

def dna_to_aa(sequence, try_frames=False):
    """
    sequence (string) DNA sequence
    will search for correct reading frame
    """
    orig_code = genetic_code(11)

    if not try_frames:
        return orig_code.translate(sequence).sequence

    sequence = DNASequence(sequence)
    translated = orig_code.translate_six_frames(sequence)
    stops = [aastring.sequence.count('*') for aastring in translated]

    return translated[stops.index(min(stops))].sequence


def all_dna_point_mutants_to_aa(wt_sequence):
    """
    wt_sequence (string) DNA sequence
    assumes starting from correct reading frame
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
            if '*' in this_sequence[:-1]:
                continue
            else:
                AA_sequences.add(this_sequence)
    return AA_sequences


def two_dna_point_mutants_to_aa(wt_sequence):
    """
    wt_sequence (string) DNA sequence
    assumes starting from correct reading frame
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
                    if '*' in this_sequence[:-1]:
                        continue
                    else:
                        AA_sequences.add(this_sequence)
    return AA_sequences

