# -*- coding: utf-8 -*-

# Analyze incidence of inhibitory or activating interactions across taxon

import settings as S
import pandas as pd
import os
import numpy as np
import pdb
import scipy.stats as st
import matplotlib.pyplot as plt
from pandas.util.testing import assert_series_equal
from collections import Counter
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

reg = S.read_cache('regulation')
reg = reg[reg['Source'] == 'BRENDA'] # don't bias with just ecocyc/excluding remainder of biocyc

ki = reg[reg['Mode'] == '-']
act = reg[reg['Mode'] == '+']

#ki = S.get_data_df('inhibiting')
#act = S.read_cache('activating')
tax = S.read_cache('TaxonomicData') # was TaxonomicData_temp

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

res = pd.DataFrame( columns = ['Type','Key','EC_number','LigandID','Compound','TotalEntries','Entropy','Summary','URL','NullEntropy','NullSummary','Literature','NumReferences','Organisms'] )

for dtype in ['ki','act']:
    merge2use = ki_merge if dtype == 'ki' else act_merge
    d2use = ki if dtype =='ki' else act
    
    for g in merge2use.groups.keys():

        if len(merge2use.groups[ g ]) >= minsize:
        
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
            res.at[ixname,'Organisms'] = ';'.join(subdf['Organism'].unique())
            
            # Also calculate the entropy of all regulators of this EC, to see if it is specific to this metabolite or related to all metabolites
            bigdf = d2use[d2use['EC_number'] == g[0]]['taxonomy'].value_counts()
            res.ix[ixname,'NullEntropy'] = norm_entropy( bigdf )
            res.ix[ixname,'NullSummary'] = summarystring( bigdf )

            

# Calculate change in normalized entropy. Careful, this is not a proper statistical measure, just a heuristic!
res['DeltaEntropy'] = res['Entropy'] - res['NullEntropy']

# Write
res.to_csv('../res/Regulation_by_taxon.csv')

# Reduce to just highly confident interactions
res_reduced = res[ res['NumReferences'] >= 10]
res_reduced = res_reduced[ res_reduced['LigandID'] != '2' ]

# Keep only data with at least 10 references and ligand id 2 (= "more")
res_reduced.to_csv('../res/Regulation_by_taxon_highconfidence.csv')

# Reduce data to only the EC's in central carbon metabolism
res_reduced_EC = res_reduced[ res_reduced['EC_number'].isin(ccm['EC']) ]
res_reduced_EC['EcoliGene'] = np.nan
for ii in res_reduced_EC.index:
    res_reduced_EC.at[ii,'EColigene'] = ccm.at[ res_reduced_EC.at[ii,'EC_number'],'EcoliGene' ]
    
# Keep only data with at least 10 references and ligand id 2 (= "more")
res_reduced_EC.to_csv('../res/Regulation_by_taxon_CCM.csv')

# Make a table indicating of the number of interactions for each species
ki_species = ki.groupby(['Organism'])
act_species = act.groupby(['Organism'])
uqspecies = set(ki_species.groups.keys()).union(act_species.groups.keys())
species_df = pd.DataFrame( columns = ['Inhibition','Activation'] )
for species in uqspecies:

    if species in ki_species.groups.keys():
        species_df.at[ species,'Inhibition'] = ki_species.groups[species].size
    else:
        species_df.at[ species,'Inhibition'] = 0
    
    if species in act_species.groups.keys():
        species_df.at[ species,'Activation'] = act_species.groups[species].size
    else:
        species_df.at[ species,'Activation'] = 0
species_df['Total'] = species_df['Activation'] + species_df['Inhibition']
species_df.to_csv('../res/Regulation_by_taxon_speciescounts.csv')


##########################################################################################
# A last little piece, compare to a prior version of the result to see if anything changed
oldres = pd.read_csv('../oldres/March2017/Regulation_by_taxon.csv',header = 0,index_col = 0)

# Compare indexes
print('New result dataframe shape:')
print(res.shape)
print('Old result dataframe shape:')
print(oldres.shape)
diffix1 = oldres.index.difference(res.index)
diffix2 = res.index.difference(oldres.index)
ixix = res.index.intersection(oldres.index)

col2use=['Type','Key','EC_number','LigandID','Compound','TotalEntries','Entropy	Summary']
restrim = res.ix[ixix,col2use].astype(str)
oldrestrim = oldres.ix[ixix,col2use].astype(str)

if restrim.equals(oldrestrim):
    print('Dataframes are equal except for the new indices.')
else:
    print('Even for common indices, the old and new results have some differences.')
    # It's likely that the entropy columns have changed slightly. Don't worry too much about this.
    diffrows = np.where(restrim!=oldrestrim)[0]
    diffcols = np.where(restrim!=oldrestrim)[1]
    colcount = dict( Counter(diffcols) )
    colcount2 = {restrim.columns[item]:colcount[item] for item in colcount.keys()}
    print('Here are the different columns:')
    print(colcount2)
    

## Also check the CCM data
#oldccm = pd.read_csv('../oldres/March2017/Regulation_by_taxon_CCM.csv',header = 0,index_col = 0)
#
#newccm = res_reduced_EC.copy().astype(str)
#newccm = newccm.drop( 'Organisms',axis = 1)
#oldccm = oldccm.astype(str)
#if res_reduced_EC.equals(oldccm):
#    print('Central carbon metabolism data is unchanged.')
#else:
#    print('There are some differences in the central carbon metabolism data.')
#    
#    # It's likely that the entropy columns have changed slightly. Don't worry too much about this.
#    diffcols = np.where(newccm!=oldccm)[1]
#    colcount = dict( Counter(diffcols) )
#    colcount2 = {newccm.columns[item]:colcount[item] for item in colcount.keys()}
#    print('Here are the different columns:')
#    print(colcount2)