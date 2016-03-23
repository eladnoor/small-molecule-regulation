#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:26:52 2016

@author: noore
"""

import pandas as pd
import settings
import os

#%%
kegg_df = pd.DataFrame.from_csv(settings.KEGG2CHEBI_FNAME,
                                header=0, index_col=None)

# remove compounds that have no ChEBI:
kegg_df = kegg_df[~kegg_df['ChEBI'].isnull()]
kegg_df.set_index('KEGG_ID', inplace=True)

# split compounds with more than one ChEBI to several rows:
tmp_chebi_df = pd.DataFrame(index=kegg_df.index, 
                            data=kegg_df['ChEBI'].apply(str.split).tolist())
kegg_df = kegg_df.join(tmp_chebi_df)
kegg_df = kegg_df.drop('ChEBI', axis=1)
kegg_df['KEGG_ID'] = kegg_df.index
kegg_df = pd.melt(kegg_df, id_vars=['KEGG_ID', 'name'], value_name='ChEBI')
kegg_df = kegg_df.drop('variable', axis=1)
kegg_df = kegg_df[~kegg_df['ChEBI'].isnull()]
kegg_df.sort_values('KEGG_ID', inplace=True)


#%%

ligand_df = settings.get_data_df('ligand_ids')
ligand_df = ligand_df[~ligand_df['chebiID'].isnull()]
ligand_df['ChEBI'] = ligand_df['chebiID'].apply(lambda s: s[6:])
ligand_df = pd.merge(ligand_df[['LigandID', 'ChEBI']],
                     kegg_df[['ChEBI', 'KEGG_ID']], how='inner', on='ChEBI')
                     
# in order not to duplicate BRENDA entries, we must have exactly one KEGG_ID
# for each LigandID. Therefore, we group by the LigandIDs and keep only the
# minimal KEGG_ID
ligand_df.drop_duplicates(subset=['LigandID'], keep='first', inplace=True)
#%%

brenda_input = [{'fname': 'ecoli_turnover', 'value_col': 'Turnover_Number'},
                {'fname': 'ecoli_ki', 'value_col': 'KI_Value'},
                {'fname': 'ecoli_km', 'value_col': 'KM_Value'},
                {'fname': 'ecoli_activating_compounds', 'value_col': None}]
                
for d in brenda_input:
    df = settings.get_data_df(d['fname'])
    if d['value_col'] is not None:    
        df.rename(columns={d['value_col']: 'Value'}, inplace=True)
        df = pd.merge(df, ligand_df, how='inner', on='LigandID')
        df = df[['EC_number', 'Value', 'KEGG_ID', 'Commentary']]
    else:
        df = pd.merge(df, ligand_df, how='inner', on='LigandID')
        df = df[['EC_number', 'KEGG_ID', 'Commentary']]
        
    df.to_csv(os.path.join(settings.CACHE_DIR, d['fname'] + '_kegg.csv'))

