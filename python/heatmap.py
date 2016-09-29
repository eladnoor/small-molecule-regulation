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
sns.set('paper', style='white')

organism = 'Escherichia coli'
SAT_FORMULA_S = r'$[S]/\left([S] + K_S\right)$'
SAT_FORMULA_M = r'$[S]/\left([S] + K_M\right)$'
SAT_FORMULA_I = r'$[S]/\left([S] + K_I\right)$'

def get_kinetic_param(name, value_col, organism='Escherichia coli'):
    k = S.read_cache(name)
    
    # filter by organsim
    k = k[k['Organism'] == organism]
    
    # remove values with unmatched ligand
    k = k[pd.notnull(k['bigg.metabolite'])]

    # remove entries lacking quantitative data
    k = k[k[value_col] > 0]

    return k[['EC_number', 'bigg.metabolite', value_col]]

def calc_sat(k, value_col, conc_df):
    # choose the minimum value among all repeats    
    k = k.groupby(['EC_number', 'bigg.metabolite'])[value_col].min().reset_index()
    
    # join data with measured concentrations
    k = k.join(conc_df, on='bigg.metabolite', how='inner')
    
    # melt table so each line will be a combination of EC, substrate/inhibitor and
    # growth condition
    k = pd.melt(k, id_vars=('EC_number', 'bigg.metabolite', value_col),
                var_name='growth condition', value_name='concentration')
    
    k['saturation'] = k['concentration'] / (k['concentration'] + k[value_col])
    k['met'] = k['bigg.metabolite'].apply(lambda x: x[:-2])
    k['met:EC'] = k['met'].str.cat(k['EC_number'], sep=':')

    return k

def calc_median_sat(k):
    """
        calculates the [S]/K_S for all matching EC-metabolite pairs,
        in log2-fold-change.
        
        Input:
            K_df    - a DataFrame with three columns: EC_number, bigg.metabolite, Value
            conc_df - a DataFrame with 
    """
    fc_med = k.groupby(('met', 'growth condition')).median()[['saturation']].reset_index()
    fc_med = fc_med.pivot('met', 'growth condition', 'saturation')
    return fc_med.sort_index(axis=0)
    
#%%
_df = pd.DataFrame.from_csv(S.ECOLI_METAB_FNAME)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

colmap = dict(map(lambda x: (x, x[:-7]), met_conc_mean.columns))
met_conc_mean.rename(columns=colmap, inplace=True)

km_raw = get_kinetic_param('km', 'KM_Value')
km = calc_sat(km_raw, 'KM_Value', met_conc_mean)

ki_raw = get_kinetic_param('ki', 'KI_Value')
ki = calc_sat(ki_raw, 'KI_Value', met_conc_mean)

ki.to_csv(os.path.join(S.RESULT_DIR, 'ki_vs_conc.csv'))
km.to_csv(os.path.join(S.RESULT_DIR, 'km_vs_conc.csv'))

#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
km_sat_med = calc_median_sat(km)
ki_sat_med = calc_median_sat(ki)

sat_joined = km_sat_med.join(ki_sat_med, how='inner', lsuffix='_sub', rsuffix='_inh')
sat_joined = sat_joined.reindex_axis(sat_joined.mean(axis=1).sort_values(axis=0).index, axis=0)

fig, ax = plt.subplots(1, 1, figsize=(12, 10))

sns.heatmap(sat_joined, ax=ax, mask=sat_joined.isnull(),
            cbar=True, vmin=0, vmax=1, annot=True, cmap='viridis')
plt.xticks(rotation=90)
plt.yticks(rotation=0)

fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation_median.svg'))


#%%
# draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
km_pivoted = km.pivot('met:EC', 'growth condition', 'saturation')
ki_pivoted = ki.pivot('met:EC', 'growth condition', 'saturation')

km_pivoted = km_pivoted.reindex_axis(km_pivoted.mean(axis=1).sort_values(axis=0).index, axis=0)
ki_pivoted = ki_pivoted.reindex_axis(ki_pivoted.mean(axis=1).sort_values(axis=0).index, axis=0)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 30))
sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
            cbar=False, vmin=0, vmax=1, cmap='viridis')
sns.heatmap(ki_pivoted, ax=ax1, mask=ki_pivoted.isnull(),
            cbar=True, vmin=0, vmax=1, annot=True, cmap='viridis')
ax0.set_title(SAT_FORMULA_M)
ax0.set_ylabel('')
ax1.set_title(SAT_FORMULA_I)
ax1.set_ylabel('')
fig.tight_layout()
fig.savefig(os.path.join(S.RESULT_DIR, 'heatmap_saturation.svg'))

#%% Compare the CDFs of the two fold-change types (for Ki and Km)

fig, axs = plt.subplots(1, 3, figsize=(7.5, 3), sharey=True)

met_intersection = set(km['bigg.metabolite']).intersection(ki['bigg.metabolite'])
km_inter = km[km['bigg.metabolite'].isin(met_intersection)]
ki_inter = ki[ki['bigg.metabolite'].isin(met_intersection)]

ax = axs[0]
np.log2(pd.melt(met_conc_mean)['value']).hist(cumulative=True, normed=1, bins=1000,
                                              histtype='step', ax=ax, linewidth=1, color='orange')
ax.set_xlim(-7, 4)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$\log_2 [S]$ (in mM)')
ax.set_ylabel(r'Cumulative distribution')
ax.set_title('Measured metabolite conc.')
ax.legend(loc='upper left')

ax = axs[1]
np.log2(1.0/km_inter['KM_Value']).hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='substrates $(K_M)$', linewidth=1)
np.log2(1.0/ki_inter['KI_Value']).hist(cumulative=True, normed=1, bins=1000,
                                 histtype='step', ax=ax, label='inhibitors $(K_I)$', linewidth=1)
ax.set_xlim(-10, 11)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$\log_2 K_S$ (in mM)')
ax.set_title(r'Measured $K_{\rm S}$ values')
ax.legend(loc='upper left')

# compare Km and Ki for the intersection of EC numbers 

ax = axs[2]
km_inter['saturation'].hist(cumulative=True, normed=1, bins=1000,
                            histtype='step', ax=ax, label='substrates', linewidth=1)
ki_inter['saturation'].hist(cumulative=True, normed=1, bins=1000,
                            histtype='step', ax=ax, label='inhibitors', linewidth=1)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel(SAT_FORMULA_S)
ax.set_title(r'Saturation levels')
ax.legend(loc='upper left')
fig.tight_layout()

fig.savefig(os.path.join(S.RESULT_DIR, 'saturation_histogram.svg'))
fig.savefig(os.path.join(S.RESULT_DIR, 'saturation_histogram.png'), dpi=600)

#%% draw violin plots of the saturation distributions
ki['type'] = 'Inhibitor'
km['type'] = 'Substrate'
data = pd.concat([ki, km])
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.violinplot(x='type', y='saturation', data=data, ax=ax)
ax.set_ylabel(SAT_FORMULA_S)

#%%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.set_xscale('log')
ax.set_yscale('log')
ki_grouped = ki[['met:EC', 'saturation']].groupby('met:EC')
km_grouped = km[['met:EC', 'saturation']].groupby('met:EC')
merge_min = ki_grouped.min().join(km_grouped.min(), how='inner', rsuffix='_km_min', lsuffix='_ki_min')
merge_max = ki_grouped.max().join(km_grouped.max(), how='inner', rsuffix='_km_max', lsuffix='_ki_max')
merge_mean = ki_grouped.mean().join(km_grouped.mean(), how='inner', rsuffix='_km_mean', lsuffix='_ki_mean')
merged = pd.concat([merge_min, merge_max, merge_mean], axis=1)
for i in merged.index:
    ax.plot(merged.loc[i, ['saturation_km_min', 'saturation_km_max']],
            merged.loc[i, ['saturation_ki_mean', 'saturation_ki_mean']],
            '-', color=(0.5, 0.5, 0.9))
    ax.plot(merged.loc[i, ['saturation_km_mean', 'saturation_km_mean']],
            merged.loc[i, ['saturation_ki_min', 'saturation_ki_max']],
            '-', color=(0.5, 0.5, 0.9))
    ax.text(merged.loc[i, ['saturation_km_mean']],
            merged.loc[i, ['saturation_ki_mean']],
            i)
ax.set_xlim(1e-4, 1e4)
ax.set_ylim(1e-4, 1e4)
ax.plot([1e-2, 1e4], [1e-2, 1e4], 'k--')
ax.set_xlabel(r'$\log_2 \left( \frac{[S]}{K_M} \right)$')
ax.set_ylabel(r'$\log_2 \left( \frac{[S]}{K_I} \right)$')
