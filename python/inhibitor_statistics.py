# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 19:12:22 2016

@author: noore
"""

import settings as S
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import bigg
import numpy as np
from matplotlib_venn import venn3

def venn3_sets(set_a, set_b, set_c, set_labels, ax):
    # order of values for Venn diagram: (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    Abc = len(set_a.difference(set_b.union(set_c)))
    aBc = len(set_b.difference(set_a.union(set_c)))
    abC = len(set_c.difference(set_a.union(set_b)))
    ABc = len(set_a.intersection(set_b).difference(set_c))
    AbC = len(set_a.intersection(set_c).difference(set_b))
    aBC = len(set_b.intersection(set_c).difference(set_a))
    ABC = len(set_a.intersection(set_b).intersection(set_c))
    venn3(subsets = (Abc, aBc, ABc, abC, AbC, aBC, ABC),
          set_labels=set_labels, ax=ax)

ki = S.read_cache('ki')
ki = ki[ki['Organism'] == 'Escherichia coli']
ki_unique = ki.groupby(['EC_number', 'bigg.metabolite']).first().reset_index()

act = S.read_cache('activating')
act = act[act['Organism'] == 'Escherichia coli']
act_unique = act.groupby(['EC_number', 'bigg.metabolite']).first().reset_index()

interactions = pd.concat([ki_unique[['EC_number', 'bigg.metabolite']],
                          act_unique[['EC_number', 'bigg.metabolite']]])

int_count_EC = interactions.groupby('EC_number').count()
int_count_EC.sort_values('bigg.metabolite', inplace=True, ascending=False)
int_count_EC.rename(columns={'bigg.metabolite': 'count(metabolites)'}, inplace=True)
int_count_EC.to_csv(os.path.join(S.RESULT_DIR, 'count_interactions_per_EC_number.csv'))

int_count_met = interactions.groupby('bigg.metabolite').count()
int_count_met.sort_values('EC_number', inplace=True, ascending=False)
int_count_met.rename(columns={'EC_number': 'count(EC)'}, inplace=True)
int_count_met.to_csv(os.path.join(S.RESULT_DIR, 'count_interactions_per_metabolite.csv'))

#%%
# compare the shadow prices of the known interaction pairs the global
# distribution

shadow_df = S.read_cache('shadow_prices_0.1')
shadow_df.columns = map(str.lower, shadow_df.columns)

# map the EC numbers to the bigg.reaction IDs
model_reactions = S.get_reaction_table_from_xls()
bigg2ec = bigg.get_reaction_df()
chebi2bigg = bigg.get_metabolite_df()

# change all reaction IDs to lower-case (apparently the case standards have changed
# since the model was published).

ki_bigg = ki_unique.join(bigg2ec, on='EC_number', how='inner')
act_bigg = act_unique.join(bigg2ec, on='EC_number', how='inner')

all_mets = set(chebi2bigg['bigg.metabolite'].unique())
inh_mets = set(ki_bigg['bigg.metabolite'].unique())
act_mets = set(act_bigg['bigg.metabolite'].unique())

all_rxns = set(bigg2ec['bigg.reaction'].unique())
inh_rxns = set(ki_bigg['bigg.reaction'].unique())
act_rxns = set(act_bigg['bigg.reaction'].unique())

regulating_mets = inh_mets.union(act_mets)
regulated_rxns = inh_rxns.union(act_rxns)

#%%
df = shadow_df.loc[regulating_mets, regulated_rxns]
df['bigg.metabolite'] = df.index # add the bigg IDs as another column, in order to "melt" the table later

df = pd.melt(df,
             id_vars=['bigg.metabolite'],
             var_name='bigg.reaction',
             value_name='shadow_price')
df = df[pd.notnull(df['shadow_price'])]
ki_sp = pd.merge(df, ki_bigg, on=['bigg.metabolite', 'bigg.reaction'], how='left')
act_sp = pd.merge(df, act_bigg, on=['bigg.metabolite', 'bigg.reaction'], how='left')

ki_minus = ki_sp.loc[pd.isnull(ki_sp['EC_number']), 'shadow_price']
ki_plus = ki_sp.loc[pd.notnull(ki_sp['EC_number']), 'shadow_price']

act_minus = act_sp.loc[pd.isnull(act_sp['EC_number']), 'shadow_price']
act_plus = act_sp.loc[pd.notnull(act_sp['EC_number']), 'shadow_price']

#%%
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].hist(ki_minus.values, bins=1000, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[0].hist(ki_plus.values, bins=200, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[0].set_xlim(-100, 0)
axs[0].set_ylim(0, 1)
axs[0].legend(['non-interacting', 'interacting'], loc='upper left')
axs[0].set_title('Inhibiting')
axs[0].set_xlabel('Shadow Price')
axs[0].set_ylabel('Cumulative Distribution')

axs[1].hist(act_minus.values, bins=10000, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[1].hist(act_plus.values, bins=267, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[1].set_xlim(-100, 0)
axs[1].set_ylim(0, 1)
axs[1].legend(['non-interacting', 'interacting'], loc='upper left')
axs[1].set_title('Activating')
axs[1].set_xlabel('Shadow Price')
axs[1].set_ylabel('Cumulative Distribution')

fig.savefig(os.path.join(S.RESULT_DIR, 'shadow_price_cdf.svg'))

#%% 
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

venn3_sets(inh_mets, act_mets, all_mets, set_labels=('Inhibitors', 'Activators', 'All metabolites'), ax=axs[0])
venn3_sets(inh_rxns, act_rxns, all_rxns, set_labels=('Inhibited', 'Activated', 'All reactions'), ax=axs[1])
fig.savefig(os.path.join(S.RESULT_DIR, 'venn.svg'))
