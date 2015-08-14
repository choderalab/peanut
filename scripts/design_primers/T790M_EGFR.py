########################################
########################################
import os.path
from primer_design import make_single_mutant
#from primer_design import make_mutant

filename = "../../nucleotide_sequences/EGFR.txt"

wt_residue = "T"
residue_number = 790
mut_residue = "M"
first_residue = 712

#######################################

with open(filename, 'r') as fi:
    wt_sequence = fi.readline()
wt_sequence = wt_sequence[:-1]

forward_primer, reverse_primer = make_single_mutant(wt_sequence, wt_residue, residue_number, mut_residue, first_res=first_residue)

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

#forward_primer, reverse_primer = make_mutant(wt_sequence, wt_residue, residue_number, mut_residue, first_res=first_residue)

#print("Forward Primer")
#print(forward_primer)
#print("Reverse Primer")
#print(reverse_primer)

