# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 09:35:06 2015

@author: dan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv 


length = pd.DataFrame.from_csv('data/all_genes_length.csv').amino_acids
abundance = pd.DataFrame.from_csv('data/Heinmann_proteomics.csv')
relative_ab = abundance.div(abundance.sum()).glucose_heinmann
#relative_ab[length.index-abundance.index] = 0
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True)

ax1.set_axis_bgcolor('#FFE6C0')
ax2.set_axis_bgcolor('#FFE6C0')

wlength = length * relative_ab
ax1.hist(length, histtype='step', color='r', lw=2, bins=80, normed=True)
ax1.axvline(length.mean(), c='b', ls=':', lw=3)
ax1.text(800, 0.00225, 'genomic', size=15, weight='bold')
#
ax2.hist(length[relative_ab.index], histtype='step', color='r', lw=2,  bins=50, 
                                             weights=relative_ab, normed=True)
ax2.axvline(wlength.sum(), c='b', ls=':', lw=3)
ax2.text(800, 0.00225, 'weighted', size=15, weight='bold')
ax1.set_ylabel('fraction of set', size=12)
ax2.set_ylabel('fraction of set', size=12)
ax2.set_xlabel('protein length (amino acids)', size=15)
ax1.set_title('protein length distribution in E. coli, complete genome (N=4303)', size=10, weight='bold')
ax1.set_xlim(0,1250)
ax1.set_ylim(0,0.0045)
ax1.set_yticks(np.arange(0, 0.005, 0.001))
ax1.set_xticks(np.arange(0, 1001, 250))
ax1.grid(color='w', ls='-', lw=1, zorder=1)
ax2.grid(color='w', ls='-', lw=1, zorder=1)
plt.tight_layout()
fig.savefig('length_hists.pdf')
