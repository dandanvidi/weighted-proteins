import pandas as pd
import csv
import numpy as np

geiger = pd.DataFrame.from_csv('data/geiger_proteomics.csv')
names = list({l.split(' ')[1][:-2]+'_geiger' for l in geiger.columns})
genes = pd.DataFrame.from_csv('data/human_length.csv').index
geiger_out = pd.DataFrame(index=genes, columns=names)

for i in np.arange(0, 33 ,3):
    cell = geiger.T[i:i+3]
    cell = 2**cell.T
    sample = cell.columns[0].split(' ')[1][:-2]+'_geiger'
    out = csv.writer(open('data/human_%s_abundance.csv'%sample, 'w'), delimiter=',')
    out.writerow(['gene', 'abundance'])

    for k, v in cell.iterrows():
        k = k.split(';')[0]        
        v = v.dropna().mean()
        if k in genes:
            if not np.isnan(v):
                geiger_out.loc[k,sample] = v
                out.writerow([k, v])


geiger_out.dropna(how='all', inplace=True)
geiger_out.to_csv('data/geiger_proteomics_.csv')