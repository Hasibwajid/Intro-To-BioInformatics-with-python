# SEQIO OBJECT
from Bio import SeqIO


#READ/PARSE FROM FASTA FILE

f= SeqIO.parse("NC_005816.fasta","fasta")
for i in f:
    #print(i.seq)
    #print(i.id)
    #print(i.description)




from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#WORKING WITH SEQ AND SEQ RECORD

str1 = Seq('AGAAGATACC')
str2 = SeqRecord(str1)
str2.id = 'id is 1'
#another way
#str2 = SeqRecord(str1,id = 'id is 1', name='self generated sequence',description = 'i am description')
str2.annotations["evidence"] = "this is a hand made seq."
str2.letter_annotations["phred_quality"] = [40, 40, 33, 44, 656, 787, 4324, 121, 112, 10]
print(str2.annotations)
print(str2.letter_annotations)







from Bio import SeqFeature

#FEATURE LOCATION in SeqFEATURE Object.

start_pos = SeqFeature.AfterPosition(5)
end_pos = SeqFeature.BetweenPosition(9, left = 6, right = 9)
myLocation = SeqFeature.FeatureLocation(start_pos,end_pos)
exect_location = SeqFeature.FeatureLocation(5,9)
print("MY location" ,myLocation)
print("My start location: ", int(myLocation.start))
print("My end Location: ", int(myLocation.end))
print(exect_location)







from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation


#ACCESS AND DO SOME OPERATION ON SEQ USING FEATURELOCATION AND SEQFEATURE

example_parent = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGA \
CTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
str_pos = int(input("enter start position: "))
end_pos = int(input("enter end position: "))
example_feature = SeqFeature(FeatureLocation(str_pos,end_pos))
feature_seq = example_parent[example_feature.location.start :  example_feature.location.end]
feature_seq_C = example_parent[example_feature.location.start :  example_feature.location.end].complement()
feature_seq_RC = example_parent[example_feature.location.start :  example_feature.location.end].reverse_complement()
feature_seq_T = example_parent[example_feature.location.start :  example_feature.location.end].translate("Standard")

print("seq of you choice is: ",feature_seq)
print("seq of you choice after complement is: ",feature_seq_C)
print("seq of you choice after reverse complement is: ",feature_seq_RC)
print("seq of you choice after Translation is: ",feature_seq_T)







#comapre two seq records via there id's allow but directly compare two seq records is not allow

record1 = SeqRecord(Seq("ACGT"), id="test")
record2 = SeqRecord(Seq("ACGT"), id="test")
#print(record1 == record2)#this operation is not allowed
print(record1.id == record2.id)
print(record1.seq == record2.seq)








from Bio.Alphabet import generic_protein
#Bio.Alphabet (to specify type of seq and format in fasta or genbank)

record = SeqRecord(Seq("ACACGAGTATGACAGGAATGATAGTAGTGACACACCCCGAGCGAGCAGTAGTAACACAGCAGCGGAGGTGTTGATGCA",generic_protein),
                       id = "2221", description = "hand made")
print(record.format("genbank"))


'''ANOTHER METHOD '''



#take a sequence input from user and covert that into FASTA format

seq = Seq(input("ENTER DNA SEQURNCE: ").upper())
myRecord = SeqRecord(seq)
myDescription = input("Enter Description for your seq: ")
myRecord.description = myDescription
myId = input("Enter id for your sequence: ")
myRecord.id = myId
print(myRecord.format('fasta'))







#take a sequence input from user and covert that into GENBANK format

seq = Seq(input("ENTER DNA/PROTEIN SEQURNCE: ").upper(),generic_protein)
#seq.alphabet = generic_protein
myRecord = SeqRecord(seq)
myId = input("Enter id for your sequence: ")
myRecord.id = myId
myDescription = input("Enter Description for your seq: ")
myRecord.description = myDescription
print(myRecord.format('genbank'))





#Parse seq from file using SeqIO object

from Bio import SeqIO
for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))






#one by one getting record using next() function
from Bio import SeqIO
record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
first_record = next(record_iterator)
print(first_record)
print(first_record.description)
second_record = next(record_iterator)
print(second_record.id)
print(second_record.description)



#parse record in file into list and then access the perticular record using list slicing technique

from Bio import SeqIO
records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))
print("Found %i records" % len(records))
print("The last record")
last_record = records[-1] #using Python's list tricks
print("The last record id is: ",last_record.id)
print("The last record sequence is: ",repr(last_record.seq))
print("The last record length is: ",len(last_record))
print("\nThe first record")
first_record = records[0] #remember, Python counts from zero
print("The first record id is: ",first_record.id)
print("The first record sequence is: ",repr(first_record.seq))
print("The first record length is: ",len(first_record))



#getting length of all the record in file and sum all of them and get total length

from Bio import SeqIO
with open("ls_orchid.gbk") as handle:
    print(sum(len(r) for r in SeqIO.parse(handle, "gb")))

#OR
with open ("ls_orchid.gbk") as handle:
    print(sum(len(r) for r in SeqIO.parse(handle,"gb")))    



#getting length of all the record in ZIP file and sum all of them and get total length

with gzip.open("ls_orchid.gbk.gz","rt") as handle:#we can write ,read from zip file using (gzip.open(filename, format)
    print(sum(len(r) for r in SeqIO.parse(handle,"gb")))






from Bio import Entrez
#FETCHING SEQURNCES FROM ENTREZ USING(Entrez.efech())
Entrez.email ="flusdhd@gmail.com"
with Entrez.efetch(db="nucleotide",rettype="fasta", retmode="text",id="6273294") as handle:
    seq_record = SeqIO.read(handle,"fasta")
    print("%s with %i "%(seq_record.id,len(seq_record.features)))



#SAME THING FOR GENBANK FORMAT
#JUST CAHNGE RETTYPE to 'gb' from 'fasta'
Entrez.email ="flusdhd@gmail.com"
with Entrez.efetch(db="nucleotide",rettype="gb", retmode="text",id="6273294") as handle:
     for seq_record in SeqIO.parse(handle, "gb"):
            print("%s %s..." % (seq_record.id, seq_record.description[:50]))
            print("Sequence length %i, %i features, from: %s"%\
                  (len(seq_record), len(seq_record.features),\
                   seq_record.annotations["source"]))










####                               #####
####                               #####
####PENDING PENDING PENDING PENDING#####
####                               #####
####                               #####
from Bio import ExPASy
from Bio import SeqIO
handle = ExPASy.get_sprot_raw("N33_HUMAN")
with ExPASy.get_sprot_raw("O23729") as handle:
    seq_record = SeqIO.parse(handle, "swiss")
    for record in seq_record:
        print(record)
        #print(seq_record.name)
        #print(seq_record.description)
        #print(repr(seq_record.seq))
        #print("Length %i" % len(seq_record))
        #print(seq_record.annotations["keywords"])



#WRITE SEQURNCES IN FASTA AND GENBANK FILE
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
rec1 = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFR"\
+"GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" \
+"NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM" \
+"SSAC", generic_protein),
id="gi|14150838|gb|AAK54648.1|AF376133_1",
description="chalcone synthase [Cucumis sativus]")
rec2 = SeqRecord(Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILK"\
+"DMVVVEIPKLGKEAAVKAIKEWGQ", generic_protein),
id="gi|13919613|gb|AAK33142.1|",
description="chalcone synthase [Fragaria vesca subsp. bracteata]")
rec3 = SeqRecord(Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFR"\
+"EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP" \
+"KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN" \
+"NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV" \
+"SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW" \
+"IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT" \
+"TGEGLEWGVLFGFGPGLTVETVVLHSVAT", generic_protein),
id="gi|13925890|gb|AAK49457.1|",
description="chalcone synthase [Nicotiana tabacum]")

my_records = [rec1, rec2, rec3]
SeqIO.write(my_records, "my_example.faa", "fasta")


#read all info of databases in ncbi using entrez search engine
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo()
result = handle.read()
handle.close()
print(result)



Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()
print(record.keys())
print(record['DbList'])
print(record.values())
print(record['DbList'])



Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo(db='protein')
record = Entrez.read(handle)
handle.close()
print('Description: ',record['DbInfo']['Description'])
print('Count: ',record['DbInfo']['Count'])
print('LastUpdate: ',record['DbInfo']['LastUpdate'])


Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.esummary(db='nlmcatalog',
                         id = '101660833')#id should be exist in db this is id of pubmed database
record = Entrez.read(handle)
handle.close()
info = record[0]['TitleMainList'][0]
print('Journal info\nid: {} \nTitle: {}'.format(record[0]['Id'],info['Title']))



