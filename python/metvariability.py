# A quick script to analyze metabolite variability as a function of regulation strength

import os, sys, numpy as np, scipy as sp, pandas as pd, pdb, matplotlib.pyplot as plt, scipy.stats as st

# Read in the interactions and metabolite concentrations
reg = pd.read_csv('../res/ecoli_interactions.csv',header = 0,index_col = 0)
conc = pd.read_csv('../data/ecoli_metabolites_kochanowski2017.csv',header = 0,index_col = 0)
conc = conc.ix[:,1:14]

reg = reg.ix[:,['bigg.metabolite','EC_number']].drop_duplicates()
mcounts = reg.groupby('bigg.metabolite').count()
mvar = conc.std(axis = 1)/conc.mean(axis = 1)

mvar.index = [item.split('_')[0] for item in mvar.index]
mcounts.index = [item.split('_')[0] if ('_d' not in item and '_t' not in item) else item for item in mcounts.index]

ixmets = set(mvar.index).intersection(mcounts.index)
df = pd.DataFrame()
df['Regulation'] = mcounts.ix[ixmets,0]
df['Variability'] = mvar.ix[ixmets]

plt.ion()
plt.figure(1)
plt.plot(mvar.ix[ixmets],mcounts.ix[ixmets,0],'o')
plt.xlabel('Temporal Variation, std/mean')
plt.ylabel('Number of Regulatory Edges')

print(st.spearmanr(mvar.ix[ixmets],mcounts.ix[ixmets,0]))