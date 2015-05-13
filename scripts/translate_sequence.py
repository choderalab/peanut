########################################
# to run this script:
# python translate_sequence.py filename_wildtype_sequence 
########################################


from convert_dna_to_aa import dna_to_aa
import sys

# read from file
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    raise IOError("python count_mutants.py filename_wildtype_sequence")

with open(filename, 'r') as fi:
    sequence = fi.readline()
sequence = sequence[:-1]

AA_sequence = dna_to_aa(sequence)

namebits = filename.split('.')
# this isn't going to work right with directories... oops
outfilename = namebits[0]+"_TRANSLATED."+namebits[1]
with open(outfilename, 'w') as fo:
    fo.write(AA_sequence)
print(AA_sequence)


