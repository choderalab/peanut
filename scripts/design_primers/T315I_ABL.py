########################################
########################################
import os.path
from peanut import primer_design

filename = "../../nucleotide_sequences/Abl.txt"

wt_residue = "T"
residue_number = 315
mut_residue = "I"
first_residue = 242

#######################################

with open(filename, 'r') as fi:
    wt_sequence = fi.readline()
wt_sequence = wt_sequence[:-1]

ABL_primer_generator = primer_design.PrimerGenerator(wt_sequence, first_res=first_residue)
forward_primer, reverse_primer = ABL_primer_generator.make_single_mutant(wt_residue, residue_number, mut_residue)

print("Forward Primer")
print(forward_primer)
print("Reverse Primer")
print(reverse_primer)
print("Length of primer: "+str(len(forward_primer))+"bp")

outfilename = "../../primers/"+wt_residue+str(residue_number)+mut_residue+"_"+filename.split('/')[-1]
if not os.path.exists(outfilename):
    with open(outfilename, 'w') as fo:
        fo.write("Forward Primer\n")
        fo.write(forward_primer)
        fo.write("\nReverse Primer\n")
        fo.write(reverse_primer)

