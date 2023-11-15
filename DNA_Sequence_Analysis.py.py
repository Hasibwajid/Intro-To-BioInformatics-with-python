from Bio.Seq import Seq
from Functions import *


# # QUESTION 1 TRANSCRIPTION AFTER CHECK IF VALID SEQUENCE
#
from Bio.Seq import Seq
dna_seq=['A','C','G','T']
dna=Seq(input('[INPUT] Enter Dna seq: ').upper())
isValid=True

for d in dna:
    if d not in dna_seq:
        isValid=False
        break
    else:
        continue

if isValid==True:
    transciption = dna.replace('T', 'U')
    print('[1] Dna to Rna: ', transciption)
else:
    print('Invalid sequence! please enter dna sequence.')



# # QUESTION 2 Translation
print('[2] CODON OF DNA ARE: ')
for i in range(0,len(transciption)-2,3):
    print('  ',transciption[i:i+3],'--',codonT[transciption[i:i+3]])

print('[3] GC CONTENT OF DNA SEQUENCE IS : ',gc_content(dna),'%')



