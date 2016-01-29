"""
Defines class PrimerGenerator which saves a DNA sequence and translates back and
forth to amino acid sequences, with methods to generate forward and reverse primers
for desired point mutations to the wild type.

Dependencies:
    Requires scikit-bio
        https://github.com/biocore/scikit-bio
"""

from skbio.sequence import genetic_code
from skbio.sequence import DNASequence

class PrimerGenerator(object):
    """
    Constructor:
    ------------
        PrimerGenerator(sequence, first_res=1)
        Arguments:
        ----------
            sequence : str
                original DNA nucleotide sequence, treated as wild type
            first_res : int
                residue id number of first residue in sequence (default = 1)

    Attributes:
    -----------
        self.sequence : string
            original DNA nucleotide sequence, input at instantiation
            treated as wild type
        self.aa_sequence : string
            wild type sequence of one-letter amino acid codes translated
            from input nucleotide sequence
        self.first_res : int
            index of the first residue in the input sequence
        self.orig_code : skbio.sequence.genetic_code
            arbitrary code object used by sci-kit bio as translator
    Methods:
    --------
        make_single_mutant()
            Determines how many nucleotide changes are required for the desired amino acid
            mutation, then constructs a primer with a minimum of 25 nucleotides, increasing 
            the length symmetrically (such that the mutant codon is centered in the primer)
            up to 45 nucleotides, using the minimum length possible to achieve acceptable
            melting temperature (78C minimum)
            Calls sub-methods:
            ------------------
                _check_melting_temp()
                    Calculates melting temp of a given primer sequence, suggests new start and
                    end indices for primer if melting temp and / or GC content is too low
                _make_mutant()
                    Finds the mutant codon, if mutation requires more than 1 nucleotide change

    """
    def __init__(self, sequence, first_res=1):
        """
        Arguments:
        ----------
            sequence : str
                original DNA nucleotide sequence, treated as wild type
            first_res : int
                residue id number of first residue in sequence (default = 1)
        """
        orig_code = genetic_code(11)
        sequence = sequence.upper()
        aa_sequence = orig_code.translate(sequence).sequence
    
        self.sequence = sequence
        self.aa_sequence = aa_sequence
        self.first_res = first_res
        self.orig_code = orig_code
        return

    def _check_melting_temp(self, primer_sequence, start_ix, end_ix, length_sequence):
        """
        Calculates melting temp of a given primer sequence, suggests new start and
        end indices for primer if melting temp and / or GC content is too low

        Called iteratively by make_single_mutant()

        Arguments:
        ----------
            primer_sequence : str
                nucleotide sequence including mutant codon to be checked
            start_ix : int
                index of the first nucleotide in the primer sequence in the full
                nucleotide sequence
            end_ix : int
                index of the final nucleotide in the primer sequence in the full
                nucleotide sequence
            length_sequence : int
                length of the full nucleotide sequence
        
        Returns:
        --------
            good_melting_temp : Bool
                True if primer meets required criteria for GC content and melting temp
            start_ix : int
                proposed index for first nucleotide in primer sequence
            end_ix : int
                proposed index for final nucleotide in primer sequence
        """
        N = len(primer_sequence)
        gc_percent = float(primer_sequence.count('g') + primer_sequence.count('c')) / N * 100.0
        mismatch_percent = 1.000 / N * 100.0
        melting_temp = 81.5 + 0.41*gc_percent - 675.0/N - mismatch_percent
        if melting_temp < 78.0 or end_ix - start_ix < 25:
            if start_ix > 0:
                start_ix +=-1
            if end_ix <length_sequence-1:
                end_ix +=1
            return False, start_ix, end_ix
        else:
            # should be actively dealing with this; for now just giving notification
            if primer_sequence[0] not in ['g','c'] and start_ix != 0:
                start_ix +=-1
                return False, start_ix, end_ix
            if primer_sequence[-1] not in ['g','c'] and end_ix != length_sequence-1:
                end_ix +=1
                return False, start_ix, end_ix
            if gc_percent < 40.0:
                print("GC out of range!")
                print(str(gc_percent)+"% GC")
            print("Melting temp: "+str(melting_temp)+"C\n")
            return True, start_ix, end_ix

    def make_single_mutant(self, wt_res,res_num,mut_res):
        """
        Determines how many nucleotide changes are required for the desired amino acid
        mutation, then constructs a primer with a minimum of 25 nucleotides, increasing 
        the length symmetrically (such that the mutant codon is centered in the primer)
        up to 45 nucleotides, using the minimum length possible to achieve acceptable
        melting temperature (78C minimum)

        DNA sequence needs to start with the first residue of the protein (no promoter, etc)
        take DNA sequence, convert to AA, define AA point mutant, find corresponding codon 
        of wt and mut, output forward and reverse primers
        DNA sequence should be only the kinase domain

        Desired mutation should require only a single nucleotide change; will print warning
        if more nucleotide changes are required

        Arguments:
        ----------
            sequence : str
                DNA sequence
            wt_res : char
                single letter amino acid code of wildtype residue to be mutated
            res_num : int 
                residue id number of residue to be mutated
            mut_res : char
                single letter amino acid code of mutant residue
        Returns:
        --------
            forward_primer : str
                nucleotide sequence
            reverse_primer : str
                nucleotide sequence
        """
        aa_sequence = self.aa_sequence
        sequence = self.sequence
        first_res = self.first_res
        orig_code = self.orig_code

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
            print("Cannot make desired mutant with a single base change")
            mut_codon = self._make_mutant(wt_codon, mut_codons)

        good_melting_temp = False
        start_ix = max(0,(res_num-first_res)*3-11)
        end_ix = min(len(sequence),(res_num+1-first_res)*3+11)

        while not good_melting_temp:
            if end_ix - start_ix > 45:
                print("Acceptable melting temp was not found")
                break
            forward_primer = sequence[start_ix:(res_num - first_res)*3]+mut_codon+sequence[(res_num+1 - first_res)*3:end_ix]
            forward_primer = forward_primer.lower()
            good_melting_temp, start_ix, end_ix = self._check_melting_temp(forward_primer, start_ix, end_ix, len(sequence))
    
        forward_sequence = DNASequence(forward_primer)
        reverse_sequence = forward_sequence.rc()
    
        reverse_primer = reverse_sequence.sequence

        return forward_primer, reverse_primer

    def _make_mutant(wt_codon, mut_codons):
        """
        Finds the mutant codon, if mutation requires more than 1 nucleotide change

        Arguments:
        ----------
            wt_codon : str
                len(wt_codon) = 3
                nucleotide codon from the wild type sequence for the residue to be mutated
            mut_codons : list(str)
                all codons that translate to desired mutant residue
        Returns:
        --------
            mut_codon : str
                codon selected from mut_codons which requires the fewest changes from
                the wild type codon
        """
        mut_codons = [DNASequence(codon) for codon in mut_codons]
        distances = [wt_codon.distance(codon) for codon in mut_codons]

        changed_bp = int(min(distances)*3)
        print("This mutant required "+str(changed_bp)+"bp modifications\n")
        # choose the codon that requires fewest changes
        return mut_codons[distances.index(min(distances))].sequence




