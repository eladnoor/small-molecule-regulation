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
model, metabolites, reactions, S = settings.get_ecoli_json()
reactions = [x.lower() for x in reactions]

S = pd.DataFrame( S, index = metabolites, columns = reactions )

Rfull = pd.DataFrame(0,index = metabolites,columns = reactions)

# Read in long-format R
R = pd.read_csv( '../cache/iJO1366_SMRN.csv',header = 0,index_col = 0 )

for ii in R.index:
    Rfull.set_value( R.at[ii,'bigg.metabolite'], R.at[ii,'bigg.reaction'], 1 )

# Histograms of metabolite and reaction regulation
gR = Rfull.sum(axis = 0)
mR =  Rfull.sum(axis = 1)

plt.figure()
plt.hist( gR,align = 'left', bins = range(np.max(gR)+1) )
plt.xlabel('Degree of Reaction')

plt.figure()
plt.hist( gR[gR!=0],align = 'left', bins = range(np.max(gR)+1) )
plt.xlabel('Degree of Reaction, Null Reactions Removed')

plt.figure()
plt.hist( mR,align = 'left', bins = range(np.max(gR)+1) )
plt.xlabel('Degree of Metabolite')

plt.figure()
plt.hist( mR[mR!=0],align = 'left', bins = range(np.max(gR)+1) )
plt.xlabel('Degree of Metabolite, Null Metabolites Removed')

# Compare degree of metabolites and degree of reactions
gS = S.abs().astype(bool).sum(axis = 0)
mS = S.abs().astype(bool).sum(axis = 1)

plt.figure()
plt.loglog(mS,mR,'o')
plt.xlabel('Degree of Metabolite,S')
plt.ylabel('Degree of Metabolite,R')

plt.figure()
plt.loglog(gS,gR,'o')
plt.xlabel('Degree of Gene,S')
plt.ylabel('Degree of Gene,R')