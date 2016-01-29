########################################
# to run this script:
# python count_mutants.py filename_wildtype_sequence 
########################################

from peanut.convert_dna_to_aa import two_dna_point_mutants_to_aa
import sys

# read from file
if len(sys.argv) == 2:
    filename = sys.argv[1]#"1_YopH_orf.txt"
else:
    raise IOError("command to run script: python count_mutants.py filename_wildtype_sequence")

with open(filename, 'r') as fi:
    wt_sequence = fi.readline()
wt_sequence = wt_sequence[:-1]

AA_sequences = two_dna_point_mutants_to_aa(wt_sequence)

#print(AA_sequences)
print(len(AA_sequences))


