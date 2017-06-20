# Analyze regulatory interactions with multiple literature references

import os, sys, numpy as np, pandas as pd, pdb, itertools
import matplotlib.pyplot as plt
plt.ion()
plt.close('all')

reg = pd.read_csv('../res/ecoli_interactions.csv',header = 0,index_col = 0)
reg = reg[reg['Source'] == 'BRENDA']

# Calculate references
reg['RefList'] = [item.split(',') if pd.notnull(item) else 0 for item in reg['Literature']]

cols = ('bigg.reaction', 'bigg.metabolite')
reglit = reg.groupby(cols)

# Make a data frame to plot the results
highc = pd.DataFrame( columns = ['NumberRefs'] )
for ii in reglit.groups.keys():

    # Concatenate all the references and take the 
    ixs = reglit.groups[ii]
    tempref = reg.ix[ixs,'RefList']
    highc.ix['<--'.join(ii),'NumberRefs'] = len(np.unique(list(itertools.chain.from_iterable(tempref))))
    
plt.hist(highc['NumberRefs'], bins = np.arange(1,highc.max().max()))