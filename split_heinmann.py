import pandas as pd
import csv

heinmann = pd.DataFrame.from_csv('data/heinmann_proteomics.csv')

for c in heinmann.columns:
    out = csv.writer(open('data/ecoli_%s_abundance.csv'%c, 'w'), delimiter=',')
    out.writerow(['gene', 'abundance'])
    for i, j in heinmann[c].iteritems():
        if j:
            out.writerow([i, j])
        