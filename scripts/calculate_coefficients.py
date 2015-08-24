########################################
# to run this script:
# python calculate_coefficients.py filename_aa_sequence 
########################################

from calculate_protein_constants import calculate_extinction_coeff, calculate_molecular_weight, calculate_absorbance
import sys

# read from file
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    raise IOError("command to run script: python calculate_coefficients.py filename_aminoacid_sequence")

with open(filename, 'r') as fi:
    aa_sequence = fi.readline()
print(aa_sequence)

secreted, cytosolic = calculate_extinction_coeff(aa_sequence)

print("")
print("Extinction Coefficient (Secreted Protein)")
print(secreted)
print("Extinction Coefficient (Cytosolic Protein)")
print(cytosolic)

molecular_weight = calculate_molecular_weight(aa_sequence)
print(" ")
print("Molecular Weight (g/mol)")
print(molecular_weight)

absorbance = calculate_absorbance(aa_sequence, secreted=secreted, molecular_weight=molecular_weight)
print(" ")
print("Absorbance")
print(absorbance)


