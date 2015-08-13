########################################
########################################

from primer_design import make_single_mutant


filename = "../nucleotide_sequences/Abl.txt"

wt_residue = "T"
residue_number = 315
mut_residue = "I"
first_residue = 242

#######################################

with open(filename, 'r') as fi:
    wt_sequence = fi.readline()
wt_sequence = wt_sequence[:-1]

forward_primer, reverse_primer = make_single_mutant(wt_sequence, wt_residue, residue_number, mut_residue, first_res=first_residue)

print("Forward Primer")
print(forward_primer)
print("Reverse Primer")
print(reverse_primer)

