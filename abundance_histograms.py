# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 09:35:06 2015

@author: dan
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv 
from scipy.stats import norm
from scipy.stats.kde import gaussian_kde
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec

class protein_length(object):
    
    def __init__(self, length, abundance, xlim=5000, ylim=1):
        
        self.xlim = xlim
        self.ylim = ylim
        
        self.length = length
        self.abundance = abundance

    def genomic_dist(self, ax, label='', draw_hist=True, draw_KDE=True, 
                     KDE_color='r', hist_color='0.6'):
        
        ax.set_axis_bgcolor('#FFE6C0')
        if draw_hist:
            n, bins, patches = ax.hist(self.length, histtype='stepfilled', 
                                       color=hist_color, edgecolor='none', lw=2, 
                                       bins=range(0, np.max(self.length) + 50, 50), 
                                       normed=True, stacked=True, alpha=0.3)

        if draw_KDE:        
            f1 = gaussian_kde(self.length)
            x1 = np.linspace(0, self.xlim, 100)
            ax.plot(x1, f1(x1), KDE_color, lw=2, label=label)

        ax.axvline(self.length.median(), c=KDE_color, ls=':', lw=3)
        ax.text(self.xlim * 0.60, self.ylim * 0.5, 'genomic\n'+
                        r'$\rm{median}=%.0f$'%self.length.median(), size=12.5, weight='bold')
        ax.set_xlim(0,self.xlim)
#        ax.set_ylim(0,self.ylim)

        ax.grid(color='w', ls='-', lw=1, zorder=0)
        ax.tick_params(color='w')
        return n, bins, patches

    def weighted_dist(self, ax, label='', draw_hist=True, draw_KDE=True,
                      KDE_color='r', hist_color='0.6'):
        
        genes = self.abundance.index & self.length.index
        weights = map(int, map(np.round, 
                               self.abundance[genes])/self.abundance[genes].mean())
               
        array = list(self.length[genes].values)
        
        warr = []
        for i, x in enumerate(weights):
            warr += [array[i]] * x
        
        ax.set_axis_bgcolor('#FFE6C0')
        
        if draw_hist:
            ax.hist(self.length[genes], histtype='stepfilled', color=hist_color, 
                     edgecolor='none', lw=2, weights=self.abundance[genes],  
                     normed=True, stacked=True, alpha=0.3,
     				 bins=range(0, np.max(warr) + 50, 50))
        if draw_KDE:
            f2 = gaussian_kde(warr)
            x2 = np.linspace(0, self.xlim, 100)
            ax.plot(x2, f2(x2), KDE_color, lw=2, label=label)

        ax.axvline(np.median(warr), c=KDE_color, ls=':', lw=3)
        ax.text(self.xlim * 0.60, self.ylim * 0.5, 'weighted\n'+r'$\rm{median}=%.0f$'%np.median(warr), size=12.5, weight='bold')        
        ax.grid(color='w', ls='-', lw=1, zorder=0)
        ax.tick_params(color='w')

def get_functional_group(group_name, extend_fname):

    j = 0    
    systematic_level = 3
    genes = []
    for row in csv.reader(open(extend_fname, 'r'), delimiter='\t'):
        if len(row) < 3:
            continue
        if j != 0:
            if row[2] != '':
                break      
            genes.append(row[systematic_level])        
  
        if row[2] == group_name:
            j = 1
    return genes




if __name__ == '__main__':
    
    data = [['ecoli','ecoli_heinmann_glucose'],
            ['yeast','yeast_nagaraj_rich'], 
            ['human','human_geiger_2012_HeLa']]

    heinmann = pd.DataFrame.from_csv('data/heinmann_proteomics.csv').columns
    data = [['ecoli', h] for h in heinmann]    

#    fig, axes = plt.subplots(2,len(data), figsize=(40,6), sharex=True, sharey=True)    
    for i, d in enumerate(data):
        fig, axes = plt.subplots(2,1, figsize=(6,8), sharex=True, sharey=True) 
        org, prot = d
        xlim=1000
        ylim=0.0035
#        extend_fname = 'data/%s/extend.csv'%o

        f_length = 'data/%s_length.csv' %org
        f_abundance = 'data/%s_%s_abundance.csv' %(org, prot)
#
        abundance = pd.DataFrame.from_csv(open(f_abundance)).abundance.dropna()
        length = pd.DataFrame.from_csv(open(f_length)).length.dropna()
#	

#
        pl = protein_length(length, abundance,
                            xlim=xlim, ylim=ylim)
#    
        n, bins, patches = pl.genomic_dist(ax=axes[0], KDE_color='r', label='genomic')
        pl.weighted_dist(ax=axes[1], KDE_color='b', label='weighted')
#
        axes[0].set_title(prot[:-9] + ' [N=%i]'%len(pl.length), weight='bold')
        axes[0].set_ylabel('fraction of set', size=15)
        axes[1].set_ylabel('fraction of set', size=15)
        axes[1].set_xlabel('protein length (AA)', size=15)
#        axes[0,0].set_yticks(np.arange(0, ylim, 0.001))

        plt.tight_layout()
        fig.savefig('res/%s_%s_hist.pdf' %(org, prot))
#
