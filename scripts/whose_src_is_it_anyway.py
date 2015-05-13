with open("../aminoacid_sequences/4_Src_orf_TRANSLATED.txt",'r') as fi:
    src4 = fi.readline()[:-1]

with open("../aminoacid_sequences/Human_SRC.txt",'r') as fi:
    human = fi.readline()[:-1]

with open("../aminoacid_sequences/Chicken_SRC.txt",'r') as fi:
    chicken = fi.readline()[:-1]

if src4.find(human) != -1:
    print("Human Src!")
elif src4.find(chicken) != -1:
    print("Chicken Src!")
else:
    print("Messed up Src?")
