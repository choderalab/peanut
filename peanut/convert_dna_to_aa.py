"""
Converts DNA nucleotide sequences to amino acid sequences
Methods can be used within a SequenceConverter instance, which saves the input
DNA sequence (can be used if multiple translations will be performed) or 
externally by directly inputing a DNA sequence

Dependencies:
    Requires scikit-bio
        https://github.com/biocore/scikit-bio
"""
from skbio.sequence import genetic_code
from skbio.sequence import DNASequence

class SequenceConverter(object):
    """
    Attributes:
        self.wt_dna_sequence : string
            DNA nucleotide sequence
            saved unmodified from sequence input to instantiation
        self.aa_sequence : string
            amino acid sequence, as found by translating self.wt_dna_sequence
    Methods:
        dna_to_aa()
            translates from the input DNA nucleotide sequence to amino acid sequence
            saves self.aa_sequence
        all_dna_point_mutants_to_aa()
            finds unique amino acid sequences that can be achieved by making a single
            nucleotide mutation and translating the new nucleotide sequence to amino
            acid sequence
            returns set of strings
        two_dna_point_mutants_to_aa()
            finds unique amino acid sequences that can be achieved by making 2 nucleotide
            mutations and translating the new nucleotide sequence to amino acid sequence
            runs using nesting loops through every nucleotide in the sequence twice : VERY SLOW
            returns set of strings
    """
    def __init__(self, wt_dna_sequence):
        """
        Arguments:
        ----------
            wt_dna_sequence : str
               DNA nucleotide sequence
        """
        self.wt_dna_sequence = wt_dna_sequence

    def dna_to_aa(self, try_frames=False):
        """
        Translates from the input DNA nucleotide sequence to amino acid sequence

        Arguments:
        ----------
            Optional:
            ---------
                try_frames : Bool
                    if True, tries 6 possible reading frames, translates all to amino
                    acids and chooses sequence with fewest stop codons
                    default = False
        Returns:
        --------
            aa_sequence : str
                sequence of one-letter amino acid codes
        """
        aa_sequence = dna_to_aa(self.wt_dna_sequence, try_frames=try_frames)
        self.aa_sequence = aa_sequence
        return aa_sequence
    
    def all_dna_point_mutants_to_aa(self):
        """
        Finds all potential sequences which can be achieved by making a single nucleotide
        mutation and translating to amino acid sequence
        Ignores mutations that lead to nonsense instead of missense mutations
        Assumes self.wt_dna_sequence starts on the correct reading frame
        
        Returns:
        --------
            aa_sequences : set of str
                each str is a unique sequence of one-letter amino acid codes
        """
        aa_sequences = all_dna_point_mutants_to_aa(self.wt_dna_sequence)
        return aa_sequences

    def two_dna_point_mutants_to_aa(self):
        """
        Finds all potential sequences which can be achieved by making 2 nucleotide
        mutations and translating to amino acid sequence
        Ignores mutations that lead to nonsense instead of missense mutations
        Assumes self.wt_dna_sequence starts on the correct reading frame
        
        Returns:
        --------
            aa_sequences : set of str
                each str is a unique sequence of one-letter amino acid codes
        """

        aa_sequences = two_dna_point_mutants_to_aa(self.wt_dna_sequence)
        return aa_sequences

def dna_to_aa(sequence, try_frames=False):
    """
    Translates from the input DNA nucleotide sequence to amino acid sequence

    Arguments:
    ----------
        sequence : str
            DNA nucleotide sequence
        Optional:
        ---------
            try_frames : Bool
                if True, tries 6 possible reading frames, translates all to amino
                acids and chooses sequence with fewest stop codons
                default = False
    Returns:
    --------
        aa_sequence : str
            sequence of one-letter amino acid codes
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
    Finds all potential sequences which can be achieved by making a single nucleotide
    mutation and translating to amino acid sequence
    Ignores mutations that lead to nonsense instead of missense mutations
    Assumes wt_sequence starts on the correct reading frame
        
    Arguments:
    ----------
        wt_sequence : str
            DNA nucleotide sequence
    Returns:
    --------
        AA_sequences : set of str
            each str is a unique sequence of one-letter amino acid codes
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
    Finds all potential sequences which can be achieved by making 2 nucleotide
    mutations and translating to amino acid sequence
    Ignores mutations that lead to nonsense instead of missense mutations
    Assumes wt_sequence starts on the correct reading frame
        
    Makes double mutants via nested loops through each nucleotide in the sequence twice,
    is therefore VERY SLOW

    Arguments:
    ----------
        wt_sequence : str
            DNA nucleotide sequence
    Returns:
    --------
        AA_sequences : set of str
            each str is a unique sequence of one-letter amino acid codes
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

