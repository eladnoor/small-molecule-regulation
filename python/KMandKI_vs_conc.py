# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 11:03:03 2016

@author: noore
"""

import os
import pandas as pd
import settings
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

_df = pd.DataFrame.from_csv(os.path.join(settings.DATA_DIR, 'ecoli_metabolites_gerosa2015.csv'),
                            index_col=0)
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

km = pd.read_csv('../cache/ecoli_km_bigg.csv', index_col = 0)
ki = pd.read_csv('../cache/ecoli_ki_bigg.csv', index_col = 0)
km.rename(columns={'Value':'K_M'}, inplace=True)
ki.rename(columns={'Value':'K_I'}, inplace=True)

km = km[km['K_M'] != -999]
ki = ki[ki['K_I'] != -999]

km_median = km.groupby('bigg.metabolite').median()
ki_median = ki.groupby('bigg.metabolite').median()

data = km_median.join(ki_median, how='outer').join(met_conc_mean.iloc[:, [0, 3]], how='outer')
#data = np.log10(data)

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
    ax.set_ylabel('mean concentration on %s [mM]' % carbon_source)
    settings.plotdiag(ax=ax, lw=0.5)

fig.savefig(os.path.join(settings.RESULT_DIR, 'KMandKI_vs_conc.svg'))