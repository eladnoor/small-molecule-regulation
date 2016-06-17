# -*- coding: utf-8 -*-

# Analyze incidence of inhibitory or activating interactions across 

import settings as S
import pandas as pd
import os
import numpy as np
import pdb
import scipy.stats as st

def norm_entropy( series ):
    # entropy, normalized by maximum value = ln(# of entries)
    norm_entropy = st.entropy( series )/series.sum()
    return norm_entropy

tax2use = 'kingdom'
minsize = 5

ki = S.read_cache('ki')
act = S.read_cache('activating')
tax = S.read_cache('TaxonomicData_temp')

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
ki = ki[pd.notnull(ki['Compound'])]
ki = ki[pd.notnull(ki['taxonomy'])]

act = act[pd.notnull(act['Compound'])]
act = act[pd.notnull(act['taxonomy'])]

# We don't want duplicate measurements of the same EC:Compound in the same organism
ki.index = [':'.join( [ki.at[row,'EC_number'],ki.at[row,'Compound'],ki.at[row,'Organism']] ) for row in ki.index]

act.index = [':'.join([act.at[row,'EC_number'], act.at[row,'Compound'], act.at[row,'Organism'] ]) for row in act.index]

ki = ki.groupby(ki.index).first()
act = act.groupby(act.index).first()

# Now do some analysis
ki_merge = ki.groupby(['EC_number','Compound'])
act_merge = act.groupby(['EC_number', 'Compound'])

res = pd.DataFrame( columns = ['Type','Key','EC_number','Compound','Entropy','FullData','URL'] )

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
            res.ix[ixname,'Compound'] = g[1]
        
            subdf = d2use.ix[ merge2use.groups[ g ],:]['taxonomy'].value_counts()
            res.ix[ixname,'Entropy'] = norm_entropy( subdf )
            res.ix[ixname,'FullData'] = ';'.join([item +':' + str(subdf.ix[item]) for item in subdf.index])
            res.ix[ixname,'URL'] = 'http://www.brenda-enzymes.org/enzyme.php?ecno=' + g[0]

# Write
res.to_csv('../res/Regulation_by_taxon.csv')