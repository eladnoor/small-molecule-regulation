#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:26:52 2016

@author: noore
"""

import pandas as pd
import settings
import os
import sys
import zipfile

bigg2chebi = pd.DataFrame.from_csv(settings.BIGG2CHEBI_FNAME)

# read the ligand ID table
ligand_df = settings.get_data_df('ligand_ids')
sys.stderr.write('Total number of Ligand IDs in BRENDA: %d\n' % ligand_df.shape[0])

# remove all compounds that have no ChEBI
ligand_df = ligand_df[~ligand_df['chebiID'].isnull()]
sys.stderr.write('Out of which have a ChEBI number: %d\n' % ligand_df.shape[0])

# remove the word "chebi:" from the beginning of the string
ligand_df['chebi'] = ligand_df['chebiID'].apply(lambda s: int(s[6:]))

# merge the ligand ID table with the bigg2chebi table, to add the BiGG
# IDs as well
ligand_df = pd.merge(ligand_df[['LigandID', 'chebi']], bigg2chebi, how='inner', on='chebi')
sys.stderr.write('Number of mappings between Ligand IDs and BiGG IDs: %d\n' % ligand_df.shape[0])
         
# in order not to duplicate BRENDA entries, we must have exactly one BiGG ID
# for each Ligand ID. Therefore, we group by the Ligand IDs and keep only the
# minimal BiGG ID
ligand_df.drop_duplicates(subset=['LigandID'], keep='first', inplace=True)
sys.stderr.write('After removing duplicates: %d\n' % ligand_df.shape[0])

brenda_input = [{'fname': 'turnover',   'value_col': 'Turnover_Number'},
                {'fname': 'ki',         'value_col': 'KI_Value'},
                {'fname': 'km',         'value_col': 'KM_Value'},
                {'fname': 'activating', 'value_col': None}]
         
brenda_zip = zipfile.ZipFile(open(settings.ECOLI_BRENDA_ZIP_FNAME, 'r'))
for d in brenda_input:
    df = pd.DataFrame.from_csv(brenda_zip.open(d['fname'] + '.csv', 'r'), header=0, index_col=None)
    if d['value_col'] is not None:    
        df.rename(columns={d['value_col']: 'Value'}, inplace=True)
        df = pd.merge(df, ligand_df, how='inner', on='LigandID')
        df = df[['EC_number', 'Value', 'bigg.metabolite', 'Commentary', 'Literature']]
    else:
        df = pd.merge(df, ligand_df, how='inner', on='LigandID')
        df = df[['EC_number', 'bigg.metabolite', 'Commentary', 'Literature']]
        
    df.to_csv(os.path.join(settings.CACHE_DIR, 'ecoli_' + d['fname'] + '_bigg.csv'))

