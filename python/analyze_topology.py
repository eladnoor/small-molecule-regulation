# Compare topology of SMRN to topology of S using pre-computing shortest-distances

import os, sys, numpy as np, pickle, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set_style('ticks')
plt.ion()
plt.close('all')

# Read in the SMRN with shortest path lengths
smrn = pd.read_csv( '../cache/iJO1366_SMRN_spl.csv',header = 0,index_col = 0 )

# Make a histogram indicating the relative distance activating/inhibiting edges travel
fig, axs = plt.subplots(2, 1, figsize=(10, 10))
pd.DataFrame.hist(smrn[smrn['Value'] == -1], column=['SPL'],
                  grid=False, bins=20, ax=axs[0], linewidth=0)
axs[0].set_xlabel('Distance Traveled in Stoichiometric Matrix, Inhibitory Edges')
axs[0].set_title('Distribution of Shortest Path Lengths in Stoichiometric Matrix')

pd.DataFrame.hist(smrn[smrn['Value'] == 1], column=['SPL'],
                  grid=False, bins=20, ax=axs[1], linewidth=0)
axs[1].set_xlabel('Distance Traveled in Stoichiometric Matrix, Activating Edges')
axs[1].set_title('Distribution of Shortest Path Lengths in Stoichiometric Matrix')

# For each metabolite, compute the "mean" distance traveled
metdist = pd.DataFrame( {'MeanSPL' : smrn.groupby('bigg.metabolite')['SPL'].mean() } )
metdist['Number'] = smrn.groupby('bigg.metabolite')['SPL'].size()

pd.DataFrame.hist( metdist[metdist['Number'] > 4], column = ['MeanSPL'], 
                    grid = False, bins = 20, linewidth = 0 )
