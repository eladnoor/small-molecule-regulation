# -*- coding: utf-8 -*-
"""
Created on Mon May 30 12:28:07 2016

@author: noore
"""

# Make a heatmap of [S]/Ki for different enzymes and conditions

import os
import pandas as pd
import numpy as np
import settings as S
import matplotlib.pyplot as plt
import seaborn as sns
sns.set('notebook', style=None)

organism = 'Escherichia coli'

def calc_fold_change(K_df, column_name, conc_df):
    """
        calculates the [S]/K_S for all matching EC-metabolite pairs,
        in log2-fold-change.
        
        Input:
            K_df    - a DataFrame with three columns: EC_number, bigg.metabolite, Value
            conc_df - a DataFrame with 
    """
    
    # remove values of -999, which simply means there is no quantitative data
    # for this Ki/Km
    k = K_df[K_df[column_name] != -999]
    
    # group by and calculate minimum over all repeats (i.e. with the same
    # reaction and the same metabolite)
    k = k.groupby(['EC_number', 'bigg.metabolite'])[column_name].min().reset_index()
    k = k.join(conc_df, on='bigg.metabolite', how='inner')
    
    # remove the _c suffix in the metabolite names
    k['met'] = k['bigg.metabolite'].apply(lambda x: x[:-2])
    k.drop('bigg.metabolite', axis=1, inplace=True)

    k.set_index(['EC_number', 'met'], inplace=True)
    k = k.div(k[column_name], axis=0)
    k.drop(column_name, axis=1, inplace=True)
    return np.log2(k)

_df = pd.DataFrame.from_csv(S.METABOLITE_CONC_FNAME)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

colmap = dict(map(lambda x: (x, x[:-7]), met_conc_mean.columns))
met_conc_mean.rename(columns=colmap, inplace=True)

km = S.read_cache('km')
ki = S.read_cache('ki')

km = km[km['Organism'] == organism]
ki = ki[ki['Organism'] == organism]

km_fc = calc_fold_change(km, 'KM_Value', met_conc_mean)
ki_fc = calc_fold_change(ki, 'KI_Value', met_conc_mean)

km_fc_med = km_fc.reset_index().groupby('met').median()
ki_fc_med = ki_fc.reset_index().groupby('met').median()
km_fc_med.sort_index(axis=0, inplace=True)
ki_fc_med.sort_index(axis=0, inplace=True)

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 10))

sns.heatmap(km_fc_med, ax=ax0, mask=km_fc_med.isnull(), cbar=False, vmin=-15, vmax=15, annot=True)
plt.xticks(rotation=90)
plt.yticks(rotation=0)

sns.heatmap(ki_fc_med, ax=ax1, mask=ki_fc_med.isnull(), cbar=True, vmin=-15, vmax=15, annot=True)
plt.xticks(rotation=90)
plt.yticks(rotation=0)

ax0.set_title('$\log_2([S]/K_M)$')
ax0.set_ylabel('')
ax1.set_title('$\log_2([S]/K_I)$')
ax1.set_ylabel('')
fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation_median.svg'))

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 30))
sns.heatmap(km_fc, ax=ax0, mask=km_fc.isnull(), cbar=False, vmin=-15, vmax=15)
sns.heatmap(ki_fc, ax=ax1, mask=ki_fc.isnull(), cbar=True, vmin=-15, vmax=15, annot=True)
ax0.set_title('$\log_2([S]/K_M)$')
ax0.set_ylabel('')
ax1.set_title('$\log_2([S]/K_I)$')
ax1.set_ylabel('')
fig.tight_layout()
fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation.svg'))

#%% Compare the CDFs of the two fold-change types (for Ki and Km)
with plt.xkcd():
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    pd.melt(km_fc)['value'].hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='$K_M$', linewidth=2)
    pd.melt(ki_fc)['value'].hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='$K_I$', linewidth=2)
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r'$\log_2 \left( \frac{[S]}{K_S} \right)$')
    ax.set_ylabel(r'Cumulative distribution')
    ax.set_title('Saturation of substrates and inhibitors')
    ax.legend(loc='upper left')
    fig.savefig(os.path.join(S.RESULT_DIR, 'saturation_histogram.svg'))
