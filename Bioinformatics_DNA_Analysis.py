from Bio.Seq import *
from Bio.SeqUtils import GC
from Bio.Data import CodonTable
#1 Count dna

dna=Seq(input("Enter: ").upper())
tempcount={'A':0,'C':0,'G':0,'T':0}
for nu in dna:
    tempcount[nu]+=1
print('COUNT(A C G T): ',tempcount['A'],tempcount['C'],tempcount['G'],tempcount['T'])


# #2 DNA TO RNA
rna=dna.replace('T','U')
print('RNA: ',rna)



# #3 reverse complement
complemetseq={'A':'T','C':'G','G':'C','T':'A'}
reverseComplement=''
for i in dna[::-1]:
    reverseComplement+=complemetseq[i]
print('REVERSE COMPLEMENT OF DNA: ',reverseComplement)



# #4 Count GC content
print('GC COUNT PERCENT(USING LOGIC): ',(dna.count('C')+dna.count('G')) / len(dna)*100,'%')

print('GC COUNT PERCENT(USING GC MTETHOD): ',GC(dna),'%')


# transcription

transc=dna.transcribe()
print('TRANSCRIPTION: ',transc)
print('TRANSCRIPTION BACK: ',transc.back_transcribe())




# Translation: convert rna or dna seq into protein sequence
print(dna.translate(table=2,to_stop=True,stop_symbol='*'))

# if to_stop=True then it with stop translation where found stop codon(TGA,UGA,TAG,TAA,UAA)

# print(transc.translate(table=2,to_stop=False,stop_symbol='***'))

for i in range(0,len(dna)-2,3):
    new=dna[i:i+3]
    print(new.translate(table='Standard'),end=' ')

print()



# CHECK START AND STOP CODON OF DIFFERENT TABLES
Vertebrate=CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
standard=CodonTable.unambiguous_dna_by_name["Standard"]
print('Start codons: ',standard.start_codons)
print('Stop codons: ',standard.stop_codons)
print(standard)
print('Start codons: ',Vertebrate.start_codons)
print('Stop codons: ',Vertebrate.stop_codons)
print(Vertebrate)
