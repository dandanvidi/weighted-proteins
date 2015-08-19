from Bio import Entrez
from Bio import SeqIO
from collections import Counter
import csv

# download genbank from NCBI and read gb file
handle=Entrez.efetch(db='nucleotide',id='U00096.3',rettype='gb')
record = SeqIO.read(handle, "gb")

# prepare csv file for amino acids distibution by ORF
out = csv.writer(open('aa_distribution_by_ORF.csv', 'w'))
AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")
out.writerow(['bnumber']+AA_LETTERS)

# count amino acids per ORF and write to csv file    
for r in record.features:
    if r.type == 'CDS':
        try:
            data = r.qualifiers
            bnumber = data['locus_tag']
            aa_dist = Counter(data['translation'][0])
            aa_dist = [aa_dist[aa] for aa in AA_LETTERS]
            out.writerow(bnumber + aa_dist)
        except KeyError:
            print "NOT FOUND: %s - %s \n" %(data['locus_tag'][0], data['note'][0]) 

