from skbio.sequence import genetic_code
from skbio.sequence import DNASequence

class PrimerGenerator(object):
    """
    Requires scikit-bio
    https://github.com/biocore/scikit-bio
    """
    def __init__(self, sequence, first_res=1):
        orig_code = genetic_code(11)
        sequence = sequence.upper()
        aa_sequence = orig_code.translate(sequence).sequence
    
        self.sequence = sequence
        self.aa_sequence = aa_sequence
        self.first_res = first_res
        self.orig_code = orig_code
        return

    def check_melting_temp(self, primer_sequence, start_ix, end_ix, length_sequence):
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
            mut_codon = self.make_mutant(wt_codon, mut_codons)

        good_melting_temp = False
        start_ix = max(0,(res_num-first_res)*3-11)
        end_ix = min(len(sequence),(res_num+1-first_res)*3+11)

        while not good_melting_temp:
            if end_ix - start_ix > 45:
                print("Acceptable melting temp was not found")
                break
            forward_primer = sequence[start_ix:(res_num - first_res)*3]+mut_codon+sequence[(res_num+1 - first_res)*3:end_ix]
            forward_primer = forward_primer.lower()
            good_melting_temp, start_ix, end_ix = self.check_melting_temp(forward_primer, start_ix, end_ix, len(sequence))
    
        forward_sequence = DNASequence(forward_primer)
        reverse_sequence = forward_sequence.rc()
    
        reverse_primer = reverse_sequence.sequence

        return forward_primer, reverse_primer

    def make_mutant(wt_codon, mut_codons):
        """
        """

        mut_codons = [DNASequence(codon) for codon in mut_codons]
        distances = [wt_codon.distance(codon) for codon in mut_codons]

        changed_bp = int(min(distances)*3)
        print("This mutant required "+str(changed_bp)+"bp modifications\n")
        # choose the codon that requires fewest changes
        return mut_codons[distances.index(min(distances))].sequence




