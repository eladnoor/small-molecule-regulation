# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:03:03 2016

@author: noore
"""

import os
import pandas as pd
import scipy
import settings as S
import matplotlib.pyplot as plt
import seaborn as sns
sns.axes_style('whitegrid')

organism = 'Escherichia coli'

_df = pd.DataFrame.from_csv(S.METABOLITE_CONC_FNAME)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

km = S.read_cache('km')
ki = S.read_cache('ki')

km = km[km['Organism'] == organism]
ki = ki[ki['Organism'] == organism]

km = km[km['KM_Value'] != -999]
ki = ki[ki['KI_Value'] != -999]

km_median = km.groupby('bigg.metabolite')['KM_Value'].median().reset_index()
ki_median = ki.groupby('bigg.metabolite')['KI_Value'].median().reset_index()

data = pd.merge(km_median, ki_median).join(met_conc_mean, on='bigg.metabolite')
data.set_index('bigg.metabolite', inplace=True)

concensus = data[~data.isnull().any(axis=1)]

fig, axs = plt.subplots(1, 2, figsize=(14, 7))
for ax, carbon_source in zip(axs, ['Glucose', 'Acetate']):
    ax.set_xscale('log')
    ax.set_yscale('log')
    xdata = concensus['%s (mean)' % carbon_source]
    ydata1 = concensus['KM_Value']
    ydata2 = concensus['KI_Value']
    for ind in concensus.index:
        ax.plot([xdata[ind], xdata[ind]], [ydata1[ind], ydata2[ind]],
                'k-', alpha=0.4, linewidth=0.5)
        ax.annotate(ind, xy=(xdata[ind], max(ydata1[ind], ydata2[ind])),
                    xycoords='data', 
                    xytext=(0, 5), textcoords='offset points',
                    ha='center')
    ax.plot(xdata, ydata1, 'r.', label='$K_M$', markersize=10)
    ax.plot(xdata, ydata2, 'b.', label='$K_I$', markersize=10)
    
    ax.legend(loc='upper left')
    ax.set_xlabel('mean concentration on %s [mM]' % carbon_source)
    ax.set_ylabel('Affinity [mM]')
    S.plotdiag(ax=ax, lw=0.5)

fig.savefig(os.path.join(S.RESULT_DIR, 'KMandKI_vs_conc.svg'))

#%%
# scatter plot comparing the K_M or K_I and abundance of each metabolite in glucose and acetate
for param in ['KM_Value', 'KI_Value']:
    fig, axs = plt.subplots(1, 2, figsize=(14, 7))
    for ax, carbon_source in zip(axs, ['Glucose', 'Acetate']):
        ax.set_xscale('log')
        ax.set_yscale('log')
        x = '%s (mean)' % carbon_source
        concensus = data[~ (data[x].isnull() | data[param].isnull())]
        xdata = concensus[x]
        ydata = concensus[param]
        for ind in concensus.index:
            ax.annotate(ind, xy=(xdata[ind], ydata[ind]),
                        xycoords='data', 
                        xytext=(0, 5), textcoords='offset points',
                        ha='center')
        concensus.plot(kind='scatter', x=x, y=param, ax=ax)
        
        ax.legend(loc='upper left')
        ax.set_xlabel('mean concentration on %s [mM]' % carbon_source)
        ax.set_ylabel('$%s$ in [mM]' % param)
        S.plotdiag(ax=ax, lw=0.5)
        r = scipy.stats.spearmanr(xdata, ydata)
        ax.annotate('$r_{spearman}$ = %.2f (p < %.3f)' % (r.correlation, r.pvalue),
                    xy=(0.04, 0.9), xycoords='axes fraction',
                    ha='left', va='top', size=15)
    
    fig.savefig(os.path.join(S.RESULT_DIR, '%s_vs_conc.svg' % param))
