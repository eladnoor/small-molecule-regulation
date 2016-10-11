# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 18:40:38 2016

@author: noore
"""
import settings
import pandas as pd
import os, re, zipfile, json, csv
ORGNAISM = 'Escherichia coli'

#with zipfile.ZipFile(open(settings.BRENDA_LIGAND_ID_FNAME, 'r')) as z:
#    brenda2chebi = pd.DataFrame.from_csv(z.open('ligand_id_table.csv'))

with zipfile.ZipFile(open(settings.BRENDA_LIGAND_ID_FNAME, 'r')) as z:
    ligand_id_table = pd.read_csv(z.open('ligand_id_table.csv'), sep=',', escapechar='\\', quotechar='"')

ligand_id_table.sort_values('LigandID', inplace=True)
# group by InChI
#missing_ligands_df = ligand_id_table.groupby('LigandID').first()
ligands_without_chebi_df = ligand_id_table[pd.isnull(ligand_id_table['chebiID'])]

# intersect this table only with Ligand IDs that are involved in some E. coli Km, Ki or act/inh

all_relevant_ligand_ids = set()
brenda_zip = zipfile.ZipFile(open(settings.BRENDA_ZIP_FNAME, 'r'))
for d in settings.BRENDA_INPUT:
    df = pd.DataFrame.from_csv(brenda_zip.open(d['fname'] + '.csv', 'r'), index_col=None)
    ligand_ids = df.loc[df['Organism'] == ORGNAISM, 'LigandID'].unique()
    all_relevant_ligand_ids.update(ligand_ids)
    
ligands_without_chebi_df = ligands_without_chebi_df[ligands_without_chebi_df['LigandID'].isin(all_relevant_ligand_ids)]

ligands_without_chebi_df.to_csv(os.path.join(settings.RESULT_DIR, 'ligands_missing_chebi_ids.csv'))
