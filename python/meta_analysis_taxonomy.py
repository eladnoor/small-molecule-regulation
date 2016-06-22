# -*- coding: utf-8 -*-

# Analyze incidence of inhibitory or activating interactions across 

import settings as S
import pandas as pd
import os
import numpy as np
import pdb
import scipy.stats as st
import matplotlib.pyplot as plt
plt.ion()
plt.close('all')

def norm_entropy( series ):
    # entropy, normalized by maximum value = ln(# of entries)
    norm_entropy = st.entropy( series )/series.shape[0]
    return norm_entropy

def summarystring( subdf ):
    # summarizes the entries in subdf
    return ';'.join([item +':' + str(subdf.ix[item]) for item in subdf.index])

tax2use = 'kingdom'
minsize = 10

ki = S.read_cache('inhibiting')
act = S.read_cache('activating')
tax = S.read_cache('TaxonomicData_temp')

# Drop entries without organism
ki = ki[pd.notnull(ki['Organism'])]
act = act[pd.notnull(act['Organism'])]

# Convert LigandID to string
ki['LigandID'] = ki['LigandID'].astype(str)
act['LigandID'] = act['LigandID'].astype(str)

# Annotate with taxonomy of choice
ki = ki[ki['Organism'].isin( tax.index )]
act = act[act['Organism'].isin( tax.index )]

ki_tax = tax.ix[ ki['Organism'], tax2use ]
ki_tax.index = ki.index
ki['taxonomy'] = ki_tax

act_tax = tax.ix[ act['Organism'], tax2use ]
act_tax.index = act.index
act['taxonomy'] = act_tax 

# Drop null values
ki = ki[pd.notnull(ki['LigandID'])]
ki = ki[pd.notnull(ki['taxonomy'])]

act = act[pd.notnull(act['LigandID'])]
act = act[pd.notnull(act['taxonomy'])]

# We don't want duplicate measurements of the same EC:LigandID in the same organism
ki.index = [':'.join( [ki.at[row,'EC_number'],ki.at[row,'LigandID'],ki.at[row,'Organism']] ) for row in ki.index]

act.index = [':'.join([act.at[row,'EC_number'], act.at[row,'LigandID'], act.at[row,'Organism'] ]) for row in act.index]

ki = ki[~ki.index.duplicated()]
act = act.groupby(act.index).first()

# Now do some analysis
ki_merge = ki.groupby(['EC_number','LigandID'])
act_merge = act.groupby(['EC_number', 'LigandID'])

res = pd.DataFrame( columns = ['Type','Key','EC_number','LigandID','Compound','TotalEntries','Entropy','Summary','URL','NullEntropy','NullSummary'] )

for dtype in ['ki','act']:
    merge2use = ki_merge if dtype == 'ki' else act_merge
    d2use = ki if dtype =='ki' else act
    
    for g in merge2use.groups.keys():

        if len(merge2use.groups[ g ]) > minsize:
        
            # Get counts for each taxon
            ixname = dtype + ':' + ':'.join(list(g))
            res.ix[ixname,'Key'] = g
            res.ix[ixname,'Type'] = dtype
            res.ix[ixname,'EC_number'] = g[0]
            res.ix[ixname,'LigandID'] = g[1]
            res.ix[ixname,'Compound'] = ';'.join(d2use.ix[ merge2use.groups[ g ],:]['Compound'].unique().astype(str))
            res.ix[ixname,'TotalEntries'] = len(merge2use.groups[ g ] )
        
            subdf = d2use.ix[ merge2use.groups[ g ],'taxonomy'].value_counts()
            res.ix[ixname,'Entropy'] = norm_entropy( subdf )
            res.ix[ixname,'Summary'] = summarystring( subdf )
            
            urladd = '#INHIBITORS' if dtype == 'ki' else '#ACTIVATING%20COMPOUND'
            res.ix[ixname,'URL'] = 'http://www.brenda-enzymes.org/enzyme.php?ecno=' + g[0] + urladd
            
            # Also calculate the entropy of all regulators of this EC, to see if it is specific to this metabolite or related to all metabolites
            bigdf = d2use[d2use['EC_number'] == g[0]]['taxonomy'].value_counts()
            res.ix[ixname,'NullEntropy'] = norm_entropy( bigdf )
            res.ix[ixname,'NullSummary'] = summarystring( bigdf )

# Calculate change in normalized entropy. Careful, this is not a proper statistical measure, just a heuristic!
res['DeltaEntropy'] = res['Entropy'] - res['NullEntropy']

# Write
res.to_csv('../res/Regulation_by_taxon.csv')

# Plot results 
plt.plot( res['NullEntropy'],res['Entropy'],'o')
plt.ylabel('Entropy')
plt.xlabel('Null Entropy')