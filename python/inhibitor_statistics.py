# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 19:12:22 2016

@author: noore
"""

import settings
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

ki = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_ki_bigg.csv'))
ki_unique = ki.groupby(['EC_number', 'bigg.metabolite']).first().reset_index()

act = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_activating_bigg.csv'))
act_unique = act.groupby(['EC_number', 'bigg.metabolite']).first().reset_index()

interactions = pd.concat([ki_unique[['EC_number', 'bigg.metabolite']],
                          act_unique[['EC_number', 'bigg.metabolite']]])

int_count_EC = interactions.groupby('EC_number').count()
int_count_EC.sort_values('bigg.metabolite', inplace=True, ascending=False)
int_count_EC.rename(columns={'bigg.metabolite': 'count(metabolites)'}, inplace=True)
int_count_EC.to_csv(os.path.join(settings.RESULT_DIR, 'count_interactions_per_EC_number.csv'))

int_count_met = interactions.groupby('bigg.metabolite').count()
int_count_met.sort_values('EC_number', inplace=True, ascending=False)
int_count_met.rename(columns={'EC_number': 'count(EC)'}, inplace=True)
int_count_met.to_csv(os.path.join(settings.RESULT_DIR, 'count_interactions_per_metabolite.csv'))

#%%
# compare the shadow prices of the known interaction pairs the global
# distribution

shadow_df = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'shadow_prices.csv'))
shadow_df.columns = map(str.lower, shadow_df.columns)

# map the EC numbers to the bigg.reaction IDs
model_reactions = settings.get_reaction_table_from_xls()
bigg2ec = model_reactions.loc[:, ['Reaction Abbreviation', 'EC Number']]
bigg2ec.rename(columns={'Reaction Abbreviation': 'bigg.reaction', 'EC Number':'EC_number'},
               inplace=True)
bigg2ec = bigg2ec.loc[~bigg2ec['EC_number'].isnull()]

# change all reaction IDs to lower-case (apparently the case standards have changed
# since the model was published).

bigg2ec['bigg.reaction'] = bigg2ec['bigg.reaction'].apply(lambda s: str(s).lower())

ki_bigg = pd.merge(ki_unique, bigg2ec, on='EC_number')[['EC_number', 'bigg.metabolite', 'bigg.reaction', 'Value']]
act_bigg = pd.merge(act_unique, bigg2ec, on='EC_number')[['EC_number', 'bigg.metabolite', 'bigg.reaction']]
regulating_mets = set(ki_bigg['bigg.metabolite'].unique()).union(act_bigg['bigg.metabolite'])
regulated_rxns = set(ki_bigg['bigg.reaction'].unique()).union(act_bigg['bigg.reaction'])

#%%
df = shadow_df.loc[regulating_mets, regulated_rxns]
df['bigg.metabolite'] = df.index

df = pd.melt(df,
             id_vars=['bigg.metabolite'],
             var_name='bigg.reaction',
             value_name='shadow_price')
ki_sp = pd.merge(df, ki_bigg, on=['bigg.metabolite', 'bigg.reaction'], how='left')
act_sp = pd.merge(df, act_bigg, on=['bigg.metabolite', 'bigg.reaction'], how='left')

ki_minus = ki_sp['shadow_price'][pd.isnull(ki_sp['EC_number'])]
ki_plus = ki_sp['shadow_price'][~pd.isnull(ki_sp['EC_number'])]

act_minus = act_sp['shadow_price'][pd.isnull(ki_sp['EC_number'])]
act_plus = act_sp['shadow_price'][~pd.isnull(ki_sp['EC_number'])]

#%%
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].hist(ki_minus.values, bins=10000, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[0].hist(ki_plus.values, bins=267, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[0].set_xlim(-200, 100)
axs[0].set_ylim(0, 1)
axs[0].legend(['non-interacting', 'interacting'], loc='upper left')
axs[0].set_title('Inhibiting')
axs[0].set_xlabel('Shadow Price')
axs[0].set_ylabel('Cumulative Distribution')

axs[1].hist(act_minus.values, bins=10000, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[1].hist(act_plus.values, bins=267, cumulative=True, histtype='step', normed=True, linewidth=2)
axs[1].set_xlim(-200, 100)
axs[1].set_ylim(0, 1)
axs[1].legend(['non-interacting', 'interacting'], loc='upper left')
axs[1].set_title('Activating')
axs[1].set_xlabel('Shadow Price')
axs[1].set_ylabel('Cumulative Distribution')

fig.savefig(os.path.join(settings.RESULT_DIR, 'shadow_price_cdf.svg'))