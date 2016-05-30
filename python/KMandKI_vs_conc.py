# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:03:03 2016

@author: noore
"""

import os
import pandas as pd
import scipy
import numpy as np
import settings
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set_style('ticks')
sns.axes_style('whitegrid')

_df = pd.DataFrame.from_csv(os.path.join(settings.DATA_DIR,
                                         'ecoli_metabolites_gerosa2015.csv'),
                            index_col=0)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

km = pd.read_csv('../cache/ecoli_km_bigg.csv', index_col = 0)
ki = pd.read_csv('../cache/ecoli_ki_bigg.csv', index_col = 0)
km.rename(columns={'Value':'K_M'}, inplace=True)
ki.rename(columns={'Value':'K_I'}, inplace=True)

km = km[km['K_M'] != -999]
ki = ki[ki['K_I'] != -999]

km_median = km.groupby('bigg.metabolite')['K_M'].median().reset_index()
ki_median = ki.groupby('bigg.metabolite')['K_I'].median().reset_index()

data = pd.merge(km_median, ki_median).join(met_conc_mean, on='bigg.metabolite')
data.set_index('bigg.metabolite', inplace=True)

concensus = data[~data.isnull().any(axis=1)]

fig, axs = plt.subplots(1, 2, figsize=(14, 7))
for ax, carbon_source in zip(axs, ['Glucose', 'Acetate']):
    ax.set_xscale('log')
    ax.set_yscale('log')
    xdata = concensus['%s (mean)' % carbon_source]
    ydata1 = concensus['K_M']
    ydata2 = concensus['K_I']
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
    settings.plotdiag(ax=ax, lw=0.5)

fig.savefig(os.path.join(settings.RESULT_DIR, 'KMandKI_vs_conc.svg'))

#%%
# scatter plot comparing the K_M or K_I and abundance of each metabolite in glucose and acetate
for param in ['K_M', 'K_I']:
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
        settings.plotdiag(ax=ax, lw=0.5)
        r = scipy.stats.spearmanr(xdata, ydata)
        ax.annotate('$r_{spearman}$ = %.2f (p < %.3f)' % (r.correlation, r.pvalue),
                    xy=(0.04, 0.9), xycoords='axes fraction',
                    ha='left', va='top', size=15)
    
    fig.savefig(os.path.join(settings.RESULT_DIR, '%s_vs_conc.svg' % param))
    #

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
ki_fc = np.log2(concensus.iloc[:, 2:].multiply(1.0/concensus['K_I'], axis='index'))
km_fc = np.log2(concensus.iloc[:, 2:].multiply(1.0/concensus['K_M'], axis='index'))

colmap = dict(map(lambda x: (x, x[:-7]), ki_fc.columns))
indmap = dict(map(lambda x: (x, x[:-2]), ki_fc.index))

ki_fc.rename(index=indmap, columns=colmap, inplace=True)
km_fc.rename(index=indmap, columns=colmap, inplace=True)
fig, (ax0, ax1) = plt.subplots(1, 2, sharey=True, figsize=(20, 10))
sns.heatmap(km_fc, ax=ax0, mask=km_fc.isnull())
sns.heatmap(ki_fc, ax=ax1, mask=km_fc.isnull())
ax0.set_title('$\log_2([S]/K_M)$')
ax0.set_ylabel('')
ax1.set_title('$\log_2([S]/K_I)$')
ax1.set_ylabel('')
fig.savefig(os.path.join(settings.RESULT_DIR, 'saturation_heatmaps.svg'))