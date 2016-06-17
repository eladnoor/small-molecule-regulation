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

def get_kinetic_param(name, value_col, organism='Escherichia coli'):
    k = S.read_cache(name)
    k = k[k['Organism'] == organism]         # filter by organsim
    k = k[pd.notnull(k['bigg.metabolite'])]  # remove values with unmatched ligand
    k = k[k[value_col] > 0]                  # remove entries lacking quantitative data
    return k[['EC_number', 'bigg.metabolite', value_col]]

def calc_fc(k, value_col, conc_df):
    # choose the minimum value among all repeats    
    k = k.groupby(['EC_number', 'bigg.metabolite'])[value_col].min().reset_index()
    
    # join data with measured concentrations
    k = k.join(conc_df, on='bigg.metabolite', how='inner')
    
    # melt table so each line will be a combination of EC, substrate/inhibitor and
    # growth condition
    k = pd.melt(k, id_vars=('EC_number', 'bigg.metabolite', value_col),
                var_name='growth condition', value_name='concentration')
    
    k['log2(saturation)'] = np.log2(k['concentration'] / k[value_col])
    k['met'] = k['bigg.metabolite'].apply(lambda x: x[:-2])
    k['met:EC'] = k['met'].str.cat(k['EC_number'], sep=':')

    return k

def calc_median_fc(k):
    """
        calculates the [S]/K_S for all matching EC-metabolite pairs,
        in log2-fold-change.
        
        Input:
            K_df    - a DataFrame with three columns: EC_number, bigg.metabolite, Value
            conc_df - a DataFrame with 
    """
    fc_med = k.groupby(('met', 'growth condition')).median()[['log2(saturation)']].reset_index()
    fc_med = fc_med.pivot('met', 'growth condition', 'log2(saturation)')
    return fc_med.sort_index(axis=0)
    
_df = pd.DataFrame.from_csv(S.METABOLITE_CONC_FNAME)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

colmap = dict(map(lambda x: (x, x[:-7]), met_conc_mean.columns))
met_conc_mean.rename(columns=colmap, inplace=True)

km_raw = get_kinetic_param('km', 'KM_Value')
km = calc_fc(km_raw, 'KM_Value', met_conc_mean)

ki_raw = get_kinetic_param('ki', 'KI_Value')
ki = calc_fc(ki_raw, 'KI_Value', met_conc_mean)

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
km_fc_med = calc_median_fc(km)
ki_fc_med = calc_median_fc(ki)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 10))

sns.heatmap(km_fc_med, ax=ax0, mask=km_fc_med.isnull(),
            cbar=False, vmin=-15, vmax=15, annot=True)
plt.xticks(rotation=90)
plt.yticks(rotation=0)

sns.heatmap(ki_fc_med, ax=ax1, mask=ki_fc_med.isnull(),
            cbar=True, vmin=-15, vmax=15, annot=True)
plt.xticks(rotation=90)
plt.yticks(rotation=0)

ax0.set_title('$\log_2([S]/K_M)$')
ax0.set_ylabel('')
ax1.set_title('$\log_2([S]/K_I)$')
ax1.set_ylabel('')
fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation_median.svg'))

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
km_pivoted = km.pivot('met:EC', 'growth condition', 'log2(saturation)')
ki_pivoted = ki.pivot('met:EC', 'growth condition', 'log2(saturation)')

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 30))
sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
            cbar=False, vmin=-15, vmax=15)
sns.heatmap(ki_pivoted, ax=ax1, mask=ki_pivoted.isnull(),
            cbar=True, vmin=-15, vmax=15, annot=True)
ax0.set_title('$\log_2([S]/K_M)$')
ax0.set_ylabel('')
ax1.set_title('$\log_2([S]/K_I)$')
ax1.set_ylabel('')
fig.tight_layout()
fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation.svg'))

#%% Compare the CDFs of the two fold-change types (for Ki and Km)

fig, axs = plt.subplots(2, 3, figsize=(15, 12), sharey=True)

ax = axs[0,0]
km['log2(saturation)'].hist(cumulative=True, normed=1, bins=1000,
                            histtype='step', ax=ax, label='substrates', linewidth=2)
ki['log2(saturation)'].hist(cumulative=True, normed=1, bins=1000,
                            histtype='step', ax=ax, label='inhibitors', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$\log_2 \left( \frac{[S]}{K_S} \right)$')
ax.set_ylabel(r'Cumulative distribution')
ax.set_title('Only for EC:met pairs with $K_s$ and conc.')
ax.legend(loc='upper left')

ax = axs[0,1]
np.log2(km['concentration']).hist(cumulative=True, normed=1, bins=1000,
                                  histtype='step', ax=ax, label='substrates', linewidth=2)
np.log2(ki['concentration']).hist(cumulative=True, normed=1, bins=1000,
                         histtype='step', ax=ax, label='inhibitors', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$\log_2 [S]$ (in mM)')
ax.set_title('Only conc. that have a $K_M$ or $K_I$')
ax.legend(loc='upper left')

ax = axs[0,2]
np.log2(1.0/km['KM_Value']).hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='$K_M$', linewidth=2)
np.log2(1.0/ki['KI_Value']).hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='$K_I$', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$-\log_2 K_S$ (in mM)')
ax.set_title('Only $K_s$ that correspond to measured met.')
ax.legend(loc='upper left')

# compare Km and Ki for the intersection of EC numbers 

EC_intersection = set(km_raw['EC_number']).intersection(ki_raw['EC_number'])
km_inter = km_raw[km_raw['EC_number'].isin(EC_intersection)]
ki_inter = ki_raw[ki_raw['EC_number'].isin(EC_intersection)]

ax = axs[1,0]
np.log2(1.0/km_inter['KM_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_M$', linewidth=2)
np.log2(1.0/ki_inter['KI_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_I$', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$-\log_2 K_S$ (in mM)')
ax.set_title('Only enzymes having both $K_M$ and $K_I$')
ax.legend(loc='upper left')


met_intersection = set(km_raw['bigg.metabolite']).intersection(ki_raw['bigg.metabolite'])
km_inter = km_raw[km_raw['bigg.metabolite'].isin(met_intersection)]
ki_inter = ki_raw[ki_raw['bigg.metabolite'].isin(met_intersection)]

ax = axs[1,1]
np.log2(1.0/km_inter['KM_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_M$', linewidth=2)
np.log2(1.0/ki_inter['KI_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_I$', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$-\log_2 K_S$ (in mM)')
ax.set_title('Only metabolites having both $K_M$ and $K_I$')
ax.legend(loc='upper left')

met_intersection = set(km_raw['bigg.metabolite']).intersection(ki_raw['bigg.metabolite'])
km_inter = km_raw[km_raw['bigg.metabolite'].isin(met_intersection)]
ki_inter = ki_raw[ki_raw['bigg.metabolite'].isin(met_intersection)]

ax = axs[1,2]
met_intersection = set(met_conc_mean.index)
km_inter = km_raw[km_raw['bigg.metabolite'].isin(met_intersection)]
ki_inter = ki_raw[ki_raw['bigg.metabolite'].isin(met_intersection)]

np.log2(1.0/km_inter['KM_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_M$', linewidth=2)
np.log2(1.0/ki_inter['KI_Value']).hist(cumulative=True, normed=1, bins=1000,
                                       histtype='step', ax=ax, label='$K_I$', linewidth=2)
ax.set_xlim(-10, 10)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$-\log_2 K_S$ (in mM)')
ax.set_title('Only for metabolites that have measured conc.')
ax.legend(loc='upper left')

fig.savefig(os.path.join(S.RESULT_DIR, 'saturation_histogram.svg'))
