import pandas as pd
import csv
#sys.path.append(os.path.expanduser('~/git/across_projects'))

heinmann = pd.DataFrame.from_csv('data/heinmann_proteomics_.csv')

b_tu_upac = {row[48:54]:row[0:5] for row in open('data/all_ecoli_genes.txt', 'r')}

for c in heinmann.columns:
    out = csv.writer(open('data/ecoli_%s_heinmann_abundance.csv'%c, 'w'), delimiter=',')
    out.writerow(['gene', 'abundance'])
    for i, j in heinmann[c].iteritems():
        if j and i in b_tu_upac:
            out.writerow([b_tu_upac[i], j])
#        