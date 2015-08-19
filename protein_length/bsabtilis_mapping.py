import csv

length_file = csv.DictReader(open('cache/bsabtilis_all_genes.txt', 'r'), delimiter='\t')
length = csv.writer(open('cache/bsubtilis_length.csv', 'w'))
length.writerow(['gene', 'length'])
for i, row in enumerate(length_file):
    length.writerow([row['Locus tag'], row['Length']])


abundance_file = csv.DictReader(open('cache/224308-Spectral_counting_B.subtili_Chi_MCP_2011.txt', 'r'), delimiter='\t')
abundance = csv.writer(open('data/bsubtilis_Chi_MCP_abundance.csv', 'w'))
abundance.writerow(['gene', 'abundance'])
for i, row in enumerate(abundance_file):
    id = row['string_external_id'].split('.')[1]
    abundance.writerow([id, row['raw_spectral_count']])


#out = csv.writer(open('data/bsubtilis_Chi_MCP_abundance.csv', 'w'))
#out.writerow(['gene', 'abundance'])
#abundance = csv.reader(open('cache/bsubtili_Chi_MCP_abundance.csv', 'r'))
#abundance.next()
#i = 0
#for row in abundance:
#    ID, ab = row
#    try:
#        out.writerow([idmap[ID], ab])
#    except KeyError:
#        print "not mapped: " + ID