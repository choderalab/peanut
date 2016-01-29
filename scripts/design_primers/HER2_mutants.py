########################################
# Run this script: python HER2_mutants.py X###X
# where X###X is the desired point mutation
########################################
import os.path
from peanut import primer_design
import sys

filename = "../../nucleotide_sequences/HER2.txt"

mutant = sys.argv[1]
print("\nMutant to create: "+mutant+"\n")
wt_residue = mutant[0]
residue_number = int(mutant[1:-1])
mut_residue = mutant[-1]
first_residue = 720

#######################################

with open(filename, 'r') as fi:
    wt_sequence = fi.readline()
wt_sequence = wt_sequence[:-1]

HER2_primer_generator = primer_design.PrimerGenerator(wt_sequence, first_res=first_residue)
forward_primer, reverse_primer = HER2_primer_generator.make_single_mutant(wt_residue, residue_number, mut_residue)

print("Forward Primer")
print(forward_primer)
print("Reverse Primer")
print(reverse_primer)
print("\nLength of primer: "+str(len(forward_primer))+"bp\n")

outfilename = "../../primers/"+wt_residue+str(residue_number)+mut_residue+"_"+filename.split('/')[-1]
if not os.path.exists(outfilename):
    with open(outfilename, 'w') as fo:
        fo.write("Forward Primer\n")
        fo.write(forward_primer)
        fo.write("\nReverse Primer\n")
        fo.write(reverse_primer)
else:
    print("\nPrimer file exists; not overwritten\n")


