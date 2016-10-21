# -*- coding: utf-8 -*-

# Analyze incidence of inhibitory or activating interactions across taxon

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

def literaturestring( subdf ):
    # Summarizes literature references
    litstring = ';'.join(subdf['Literature'])
    litstring2 = ''.join(litstring.split(' '))
    uqlit = np.unique( litstring2.split(';') )
    return len(uqlit),';'.join(uqlit)

# Set some parameters
tax2use = 'kingdom'
minsize = 10

# Read in central carbon metabolism reactions
ccm = S.read_cache('CCM_Reactions')
ccm['EcoliGene'] = ccm.index
ccm.index = ccm['EC']

ki = S.read_cache('inhibiting')
act = S.read_cache('activating')
tax = S.read_cache('TaxonomicData')

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
act = act[~act.index.duplicated()]

# Remove instances of inhibition where "no inhibition" is mentioned
noinhib_ix = [item for item in ki.index if 'no inhibition' not in str(ki.at[item,'Commentary']).lower() ]

ki = ki.ix[noinhib_ix,:]

# Now do some analysis
ki_merge = ki.groupby(['EC_number','LigandID'])
act_merge = act.groupby(['EC_number', 'LigandID'])

res = pd.DataFrame( columns = ['Type','Key','EC_number','LigandID','Compound','TotalEntries','Entropy','Summary','URL','NullEntropy','NullSummary','Literature','NumReferences'] )

for dtype in ['ki','act']:
    merge2use = ki_merge if dtype == 'ki' else act_merge
    d2use = ki if dtype =='ki' else act
    
    for g in merge2use.groups.keys():

        if len(merge2use.groups[ g ]) > minsize:
        
            # Get counts for each taxon
            ixname = dtype + ':' + ':'.join(list(g))
            res.at[ixname,'Key'] = g
            res.at[ixname,'Type'] = dtype
            res.at[ixname,'EC_number'] = g[0]
            res.at[ixname,'LigandID'] = g[1]
            #res.at[ixname,'Compound'] = ';'.join(d2use.ix[ merge2use.groups[ g ],:]['Compound'].unique().astype(str))
            
            res.at[ixname,'TotalEntries'] = len(merge2use.groups[ g ] )
        
            subdf = d2use.ix[ merge2use.groups[ g ],:]
            res.at[ixname,'Compound'] = ';'.join(subdf['LigandName'].unique())
            res.at[ixname,'Entropy'] = norm_entropy( subdf['taxonomy'].value_counts() )
            res.at[ixname,'Summary'] = summarystring( subdf['taxonomy'].value_counts() )
            
            urladd = '#INHIBITORS' if dtype == 'ki' else '#ACTIVATING%20COMPOUND'
            res.at[ixname,'URL'] = 'http://www.brenda-enzymes.org/enzyme.php?ecno=' + g[0] + urladd
            
            # Get the literature references
            uqlit,uqlitstring = literaturestring( subdf )
            res.at[ixname,'Literature'] = uqlitstring
            res.at[ixname,'NumReferences'] = uqlit
            
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

# Reduce data to only the EC's in central carbon metabolism
res_reduced = res[ res['EC_number'].isin(ccm['EC']) ]
res_reduced['EcoliGene'] = np.nan
for ii in res_reduced.index:
    res_reduced.at[ii,'EColigene'] = ccm.at[ res_reduced.at[ii,'EC_number'],'EcoliGene' ]
    
# Keep only data with at least 10 references and ligand id 2 (= "more")
res_reduced = res_reduced[ res_reduced['NumReferences'] >= 10]
res_reduced = res_reduced[ res_reduced['LigandID'] != '2' ]
res_reduced.to_csv('../res/Regulation_by_taxon_CCM.csv')