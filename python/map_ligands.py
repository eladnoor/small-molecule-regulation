#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:26:52 2016

@author: noore
"""

import pandas as pd
import settings
import zipfile
import kegg
import bigg

# make a concensus table of the BRENDA ligand IDs, ChEBIs, BiGG and KEGG.

# read the ligand ID table and remove the word "chebi:" from the beginning of the string
brenda2chebi = settings.get_data_df('ligand_ids')[['LigandID', 'chebiID']]
brenda2chebi = brenda2chebi.rename(columns={'chebiID': 'ChEBI'}).set_index('LigandID')
brenda2chebi['ChEBI'] = brenda2chebi['ChEBI'].apply(lambda s: s[6:] if pd.notnull(s) else pd.np.nan)

# join the bigg IDs into the ligand table (using the ChEBIs)
chebi2bigg = bigg.get_metabolite_df()
ligand_df = brenda2chebi.join(chebi2bigg, how='left', on='ChEBI')

# combine the mapping from CIDs to ChEBIs with the mapping from ligand_IDs
# to ChEBIs to get a direct mapping from BRENDA to KEGG compound IDs
chebi2kegg = kegg.get_kegg_df()
chebi2kegg.rename(columns={'name': 'Compound'}, inplace=True)
chebi2kegg = chebi2kegg.groupby('ChEBI').first()
ligand_df = ligand_df.join(chebi2kegg, how='left', on='ChEBI')

brenda_zip = zipfile.ZipFile(open(settings.BRENDA_ZIP_FNAME, 'r'))
for d in settings.BRENDA_INPUT:
    df = pd.DataFrame.from_csv(brenda_zip.open(d['fname'] + '.csv', 'r'), index_col=None)
    df = df.join(ligand_df, how='left', on='LigandID')
    settings.write_cache(d['fname'], df)
