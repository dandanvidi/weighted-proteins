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
    
    def __init__(self, length, abundance, bins=10, xlim=5000, ylim=1):
        
        self.bins = bins
        self.xlim = xlim
        self.ylim = ylim
        
        self.length = length
        self.abundance = abundance

    def genomic_dist(self, ax, label='', draw_hist=True, draw_KDE=True, 
                     KDE_color='r', hist_color='0.6'):
        
        ax.set_axis_bgcolor('#FFE6C0')
        if draw_hist:
            ax.hist(self.length, histtype='stepfilled', color=hist_color, edgecolor='none', 
                     lw=2, bins=self.bins, normed=True, alpha=0.3)

        if draw_KDE:        
            f1 = gaussian_kde(self.length)
            x1 = np.linspace(0, self.xlim, 100)
            ax.plot(x1, f1(x1), KDE_color, lw=2, label=label)

        ax.axvline(self.length.mean(), c=KDE_color, ls=':', lw=3)
#        ax.text(self.xlim * 0.7, self.ylim * 0.5, 'genomic\n'+
#                        '  $\mu=%.0f$'%self.length.mean(), size=12.5, weight='bold')
        ax.set_xlim(0,self.xlim)
        ax.set_ylim(0,self.ylim)

        ax.grid(color='w', ls='-', lw=1, zorder=0)
        ax.tick_params(color='w')

    def weighted_dist(self, ax, label='', draw_hist=True, draw_KDE=True,
                      KDE_color='r', hist_color='0.6'):
        
        wlength = self.length * self.abundance / self.abundance.sum()
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
                     edgecolor='none', lw=2,  
                     bins=np.round(self.bins * len(genes)*1.0/len(self.length))+1, 
                     weights=self.abundance[genes], normed=True, alpha=0.3)
     
        if draw_KDE:
            f2 = gaussian_kde(warr)
            x2 = np.linspace(0, self.xlim, 100)
            ax.plot(x2, f2(x2), KDE_color, lw=2, label=label)
        ax.axvline(wlength.sum(), c=KDE_color, ls=':', lw=3)
#        ax.text(700, 0.002, 'weighted\n'+'   $\mu=%.0f$'%wlength.sum(), size=12.5, weight='bold')        
        ax.grid(color='w', ls='-', lw=1, zorder=0)
        ax.tick_params(color=KDE_color)

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

    folder = 'ecoli'
    xlim=1000
    ylim=0.007
    
    extend_fname = 'data/%s/extend.csv'%folder
    f_length = 'data/%s/length.csv'%folder
    f_abundance = 'data/%s/abundance.csv'%folder
    
    abundance = pd.DataFrame.from_csv(open(f_abundance)).abundance
    length = pd.DataFrame.from_csv(open(f_length)).length

    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True)
    
    pl = protein_length(length, abundance,
                        bins=40, xlim=xlim, ylim=ylim)
        
    pl.genomic_dist(ax=ax1, KDE_color='r', label='genomic')
    pl.weighted_dist(ax=ax2, KDE_color='b', label='weighted')
    
    genes = get_functional_group('Ribosome', extend_fname)    
    sub_ab = abundance[genes].dropna()
    sub_len = length[genes].dropna()

    sub_pl = protein_length(sub_len, sub_ab, 
                        bins=40, xlim=xlim, ylim=ylim)


    sub_pl.genomic_dist(ax=ax1, draw_hist=False, KDE_color='k', label='ribosomal')
    sub_pl.weighted_dist(ax=ax2, draw_hist=False, KDE_color='k', label='ribosomal')



    no_sub_ab = abundance[abundance.index-sub_ab.index].dropna()
    no_sub_len = length[length.index-sub_ab.index].dropna()

    no_sub_pl = protein_length(no_sub_len, no_sub_ab, 
                        bins=40, xlim=xlim, ylim=ylim)

    no_sub_pl.genomic_dist(ax=ax1, draw_hist=False, KDE_color='g', label='NO ribosomal')
    no_sub_pl.weighted_dist(ax=ax2, draw_hist=False, KDE_color='g', label='NO ribosomal')
    
    
    


    ax1.legend()
    ax2.legend()
    
    ax1.set_ylabel('fraction of set', size=12)
    ax2.set_ylabel('fraction of set', size=12)
    ax2.set_xlabel('protein length (amino acids)', size=12)
    ax1.set_yticks(np.arange(0, ylim, 0.001))
    plt.tight_layout()
    fig.savefig('res/%s/length_hists.svg'%folder)


