import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
abundance = pd.DataFrame.from_csv(open('data/ecoli_glucose_heinmann_abundance.csv')).abundance.dropna()
length = pd.DataFrame.from_csv(open('cache/ecoli_length.csv')).length.dropna()

genes = abundance.index & length[length<=600].index

x = length[genes]
y = abundance[genes] / abundance.sum()


fig = plt.figure()
ax = plt.axes()

ax.scatter(x,y, alpha=0.5, s=30, edgecolor='none', zorder=3)
a, b, r, p, std = linregress(x, np.log10(y))
plt.plot(x, 10**(a*x+b), 'r-', zorder=10)
ax.set_xlabel('length (AA)', size=15)
ax.set_ylabel('fraction of proteome', size=15)

ax.set_xlim(0,600)
#ax.set_ylim(0,0.01)

ax.grid(color='w', ls='-', lw=1, zorder=0)
ax.tick_params(color='w')
ax.set_axis_bgcolor((0.95, 0.9, 0.92))
ax.set_title('batch growth on glucose (heinmann)', size=12, weight='bold')
ax.set_yscale('log')
ax.set_ylim(1e-8,1e-1)
plt.savefig('res/length_abundance_correlation.pdf')
