########################################
# to run this script:
# python calculate_extinction_coefficient.py filename_wildtype_sequence 
########################################

from peanut.calculate_protein_constants import calculate_extinction_coeff
import sys

# read from file
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    raise IOError("command to run script: python calculate_extinction_coefficient.py filename_aminoacid_sequence")

with open(filename, 'r') as fi:
    aa_sequence = fi.readline()

secreted, cytosolic = calculate_extinction_coeff(aa_sequence)

#print(AA_sequences)
print("Extinction Coefficient (Secreted Protein)")
print(secreted)
print("Extinction Coefficient (Cytosolic Protein)")
print(cytosolic)


