import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv, os 
from scipy.stats.kde import gaussian_kde

class protein_length(object):
    '''
        Probability distributions of protein lenght across different organisms.
        Protein length is calculated in amino acids (AA), based on the 
        coding sequence in the genome.

        ABUNDANCE WEIGHTED PDF:
        
            -kernel-density estimates using Gaussian kernels
            
            -histograms with 50 AA bin width. Histograms are normalized 
            such that the integral of the histograms will sum to 1 
            to form a probability density

    '''
    
    def __init__(self, length, abundance, xlim=5000, ylim=None):
        
        self.xlim = xlim
        self.ylim = ylim
        
        self.length = length #protein length in amino acids
        self.abundance = abundance #protein abundance - realtive within samples

    def genomic_dist(self, ax, label='', draw_hist=True, draw_KDE=True, 
                     KDE_color='r', hist_color='0.6'):
        '''
            params:
                - ax: matplotlib axis to draw plot in
                - draw_hist: draws normalized histogram with 50 AA bins
                - draw_KDE: draws gaussian based kernel density estimate of data
        '''
        
        ax.set_axis_bgcolor('#FFE6C0')
        if draw_hist:
            ax.hist(self.length, histtype='stepfilled', 
                                       color=hist_color, edgecolor='none', lw=2, 
                                       bins=range(0, np.max(self.length) + 50, 50), 
                                       normed=True, stacked=True, alpha=0.3)

        if draw_KDE:        
            f1 = gaussian_kde(self.length)
            x1 = np.linspace(0, self.xlim, 100)
            ax.plot(x1, f1(x1), KDE_color, lw=2, label=label)

        # add some to plot
        ax.axvline(self.length.median(), c=KDE_color, ls=':', lw=3)
        ax.text(self.xlim * 0.60, self.ylim * 0.5, 'genomic\n'+
                        r'$\rm{median}=%.0f$'%self.length.median(), size=12.5, weight='bold')
        ax.set_xlim(0,self.xlim)
        if self.ylim:
            ax.set_ylim(0,self.ylim)
        ax.grid(color='w', ls='-', lw=1, zorder=0)
        ax.tick_params(color='w')

    def weighted_dist(self, ax, label='', draw_hist=True, draw_KDE=True,
                      KDE_color='r', hist_color='0.6'):
        '''
            params:
                - ax: matplotlib axis to draw plot in
                - draw_hist: draws normalized histogram with 50 AA bins
                - draw_KDE: draws gaussian based kernel density estimate of data
        '''
        
        genes = self.abundance.index & self.length.index

        # manually generate weigthed array in order to later plot KDE
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

        # add some to plot
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
    
    dirListing = sorted(os.listdir('data'))
    dirListing = [x for x in dirListing if not x.startswith('.')]
    
    for f in dirListing:
        org, fname = f.split('_', 1)
        fname = fname.rsplit('_', 1)[0]
        if not os.path.exists('res/%s_%s_hist.jpg' %(org, fname)):
            f_length = 'cache/%s_length.csv' %org
            fig, axes = plt.subplots(2,1, figsize=(6,8), sharex=True, sharey=True) 
            xlim=1000
            ylim=None
    #
            abundance = pd.DataFrame.from_csv(open('data/'+f)).abundance.dropna()
            length = pd.DataFrame.from_csv(open(f_length)).length.dropna()
    #	
    
    #
            pl = protein_length(length, abundance,
                                xlim=xlim, ylim=ylim)
    #    
            pl.genomic_dist(ax=axes[0], KDE_color='r', label='genomic')
            pl.weighted_dist(ax=axes[1], KDE_color='b', label='weighted')
    #   
            axes[0].set_title(fname + ' [N=%i]'%len(pl.length), weight='bold')
            axes[0].set_ylabel('fraction of set', size=15)
            axes[1].set_ylabel('fraction of set', size=15)
            axes[1].set_xlabel('protein length (AA)', size=15)
    
            plt.tight_layout()
            fig.savefig('res/%s_%s_hist.jpg' %(org, fname))
            plt.close()
#
