# Make some simple plots using the reconstructed SMRN

import os
import pandas as pd
import settings
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
plt.ion()
plt.close('all')
sns.set_style('ticks')

# Read BIGG model
model, S = settings.get_ecoli_json()

# Read in long-format R
R = pd.read_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN.csv'), header=0, index_col=0)
R_grouped = R.groupby(['bigg.reaction', 'bigg.metabolite']).sum().reset_index()
R_grouped['count'] = 1

rxn_df = pd.DataFrame(index=S.columns)
rxn_df['Degree of Reaction, R'] = R_grouped[['bigg.reaction', 'count']].groupby('bigg.reaction').sum()
rxn_df['Degree of Reaction, S'] = S.abs().astype(bool).sum(axis=0)

met_df = pd.DataFrame(index=S.index)
met_df['Degree of Metabolite, R'] = R_grouped[['bigg.metabolite', 'count']].groupby('bigg.metabolite').sum()
met_df['Degree of Metabolite, S'] = S.abs().astype(bool).sum(axis=1)

#%% Histograms of metabolite and reaction regulation
fig, axs = plt.subplots(2, 3, figsize=(14, 10))
pd.DataFrame.hist(rxn_df.fillna(0), column=['Degree of Reaction, R'],
                  grid=False, bins=20, ax=axs[0, 0], linewidth=0)
axs[0, 0].set_xlabel('Degree of Reaction')

pd.DataFrame.hist(rxn_df, column=['Degree of Reaction, R'],
                  grid=False, bins=20, ax=axs[0, 1], linewidth=0)
axs[0, 1].set_xlabel('Degree of Reaction, Null Reactions Removed')

pd.DataFrame.hist(met_df.fillna(0), column=['Degree of Metabolite, R'],
                  grid=False, bins=20, ax=axs[1, 0], linewidth=0)
axs[1, 0].set_xlabel('Degree of Metabolite')

pd.DataFrame.hist(met_df, column=['Degree of Metabolite, R'],
                  grid=False, bins=20, ax=axs[1, 1], linewidth=0)
axs[1, 1].set_xlabel('Degree of Metabolite, Null Metabolites Removed')

# Compare degree of metabolites and degree of reactions
rxn_df.plot(kind='scatter', x='Degree of Reaction, R', y='Degree of Reaction, S',
            ax=axs[0, 2], linewidth=0)
met_df.plot(kind='scatter', x='Degree of Metabolite, R', y='Degree of Metabolite, S',
            ax=axs[1, 2], linewidth=0)

fig.savefig(os.path.join(settings.RESULT_DIR, 'degree_plots.svg'))