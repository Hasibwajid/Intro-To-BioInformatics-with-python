from Bio import Entrez


# 1st
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo()
result = handle.read()
handle.close()
print(result)



#2nd
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()
print(record.keys())
print(record['DbList'])
print(record.values())
print(record['DbList'])


# 3rd
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.einfo(db='protein')
record = Entrez.read(handle)
handle.close()
print('Description: ',record['DbInfo']['Description'])
print('Count: ',record['DbInfo']['Count'])
print('LastUpdate: ',record['DbInfo']['LastUpdate'])

# 4th
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.esummary(db='nlmcatalog',
                         id = '101660833')#id should be exist in db this is id of pubmed database
record = Entrez.read(handle)
handle.close()
info = record[0]['TitleMainList'][0]
print('Journal info\nid: {} \nTitle: {}'.format(record[0]['Id'],info['Title']))


# 5th
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.esearch(db='pubmed',
                         term = 'biopython')
result = Entrez.read(handle)
handle.close()
print('Description: ',result['DbInfo']['Description'])
print('Count: ',(result['DbInfo']['Count']))
print('LastUpdate: ',(result['DbInfo']['LastUpdate']))



# 6th
'''
EFetch is what you use when you want to retrieve a full record from Entrez. This covers several possible databases '''

Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.efetch(db='nucleotide',
                       id = 'EU490707', rettype='gb',
                       retmode = 'text')
print(handle.read())
handle.close()



# 7th

# EGQuery provide counts for a search term in each of entre database(i.e a global query). 

Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.egquery(term='biopython') #here can be any database name like pubmed or any other
recode = Entrez.read(handle)
handle.close()
for row in recode['eGQueryResult']:
    print(row['DbName'],row['Count'])


# 8th
Entrez.email = 'mh9590965@gmail.com'
handle = Entrez.espell(term='biopythooon') #here can be any database name like pubmed or any other
recode = Entrez.read(handle)
recode['Query']
'biopythooon'
print(recode['CorrectedQuery'])
