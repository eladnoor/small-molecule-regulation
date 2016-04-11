#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:26:52 2016

@author: noore

This script is meant to map all the BRENDA data into E. coli b-numbers
and compound KEGG IDs, which are more standard IDs for handling such
data. Metabolite mapping to KEGG is done vie ChEBI numbers, and EC
numbers are mapped to genes using the KEGG API.
"""

import pandas as pd
import settings
import os
import kegg

kegg_df = kegg.get_kegg_df()
ec_df = kegg.get_ec_df()

# combine the mapping from CIDs to ChEBIs with the mapping from ligand_IDs
# to ChEBIs to get a direct mapping from BRENDA to KEGG compound IDs
ligand_df = settings.get_data_df('ligand_ids')
ligand_df = ligand_df[~ligand_df['chebiID'].isnull()]
ligand_df['ChEBI'] = ligand_df['chebiID'].apply(lambda s: s[6:])
ligand_df = pd.merge(ligand_df[['LigandID', 'ChEBI']],
                     kegg_df[['ChEBI', 'KEGG_ID', 'name']], how='inner', on='ChEBI')
                     
# in order not to duplicate BRENDA entries, we must have exactly one KEGG_ID
# for each LigandID. Therefore, we group by the LigandIDs and keep only the
# minimal KEGG_ID
ligand_df.drop_duplicates(subset=['LigandID'], keep='first', inplace=True)

brenda_input = [{'fname': 'ecoli_turnover', 'value_col': 'Turnover_Number'},
                {'fname': 'ecoli_ki', 'value_col': 'KI_Value'},
                {'fname': 'ecoli_km', 'value_col': 'KM_Value'},
                {'fname': 'ecoli_activating_compounds', 'value_col': None}]
                
for d in brenda_input:
    df = settings.get_data_df(d['fname'])
    df = df.join(ec_df, on='EC_number')
    df = pd.merge(df, ligand_df, how='inner', on='LigandID')
    if d['value_col'] is not None:    
        df.rename(columns={d['value_col']: 'Value'}, inplace=True)
        df = df[['Value', 'EC_number', 'b_number', 'KEGG_ID', 'name', 'Commentary']]
    else:
        df = df[['EC_number', 'b_number', 'KEGG_ID', 'name', 'Commentary']]
        
    df.to_csv(os.path.join(settings.CACHE_DIR, d['fname'] + '_kegg.csv'))

