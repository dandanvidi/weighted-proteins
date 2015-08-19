import logging
import sys
import csv
import argparse

def syntax():
    print "input error"
    sys.exit(-1)

def extend(hierarchy_fname, mapping_fname, extend_fname):
    
    mapping = {'Not mapped': 'Not Mapped:NotMapped'}
    for row in csv.reader(open(mapping_fname, 'r'), delimiter='\t'):
        systematic, gene, KO = row
        if KO in mapping:
            logging.debug('KO %s is mapped to two different genes' % KO)
        mapping[KO] = systematic

    output = csv.writer(open(extend_fname, 'w'), delimiter='\t')

    KO_level = 4

    for row in csv.reader(open(hierarchy_fname, 'r'), delimiter='\t'):
        if len(row) == KO_level:
            KO = row[-1]
            if KO not in mapping:
                logging.debug('KO %s appears in the hierarchy but is not mapped' % KO)
                continue
            row = row[:-1] + [mapping[KO]]
            
        output.writerow(row)
    
if __name__ == '__main__':

    heirary_fname = 'data/KO_gene_hierarchy_general.tms'
    mapping_fname = 'data/ecoli/eco_mapping.csv'
    extend_fname = 'data/ecoli/extend.csv'
    
    extend(heirary_fname, mapping_fname, extend_fname)
