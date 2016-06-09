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

shadow_df = S.read_cache('shadow_prices')
shadow_df.columns = map(str.lower, shadow_df.columns)

# map the EC numbers to the bigg.reaction IDs
model_reactions = S.get_reaction_table_from_xls()
bigg2ec = bigg.get_reaction_df()

# change all reaction IDs to lower-case (apparently the case standards have changed
# since the model was published).

ki_bigg = ki_unique.join(bigg2ec, on='EC_number', how='inner')
act_bigg = act_unique.join(bigg2ec, on='EC_number', how='inner')
regulating_mets = set(ki_bigg['bigg.metabolite'].unique()).union(act_bigg['bigg.metabolite'].unique())
regulated_rxns = set(ki_bigg['bigg.reaction'].unique()).union(act_bigg['bigg.reaction'].unique())

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
with plt.xkcd():
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].hist(ki_minus.values, bins=10000, cumulative=True, histtype='step', normed=True, linewidth=2)
    axs[0].hist(ki_plus.values, bins=267, cumulative=True, histtype='step', normed=True, linewidth=2)
    axs[0].set_xlim(-100, 10)
    axs[0].set_ylim(0, 1)
    axs[0].legend(['non-interacting', 'interacting'], loc='upper left')
    axs[0].set_title('Inhibiting')
    axs[0].set_xlabel('Shadow Price')
    axs[0].set_ylabel('Cumulative Distribution')
    
    axs[1].hist(act_minus.values, bins=10000, cumulative=True, histtype='step', normed=True, linewidth=2)
    axs[1].hist(act_plus.values, bins=267, cumulative=True, histtype='step', normed=True, linewidth=2)
    axs[1].set_xlim(-100, 10)
    axs[1].set_ylim(0, 1)
    axs[1].legend(['non-interacting', 'interacting'], loc='upper left')
    axs[1].set_title('Activating')
    axs[1].set_xlabel('Shadow Price')
    axs[1].set_ylabel('Cumulative Distribution')

    fig.savefig(os.path.join(S.RESULT_DIR, 'shadow_price_cdf.png'))