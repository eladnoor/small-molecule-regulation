#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:26:52 2016

@author: noore
"""

import pandas as pd
import settings
import zipfile
import os
from kegg import KEGG
from bigg import BiGG

# make a concensus table of the BRENDA ligand IDs, ChEBIs, BiGG and KEGG.

def map_brenda_to_chebi():
    # For all the Ligands that have an InChI but are not mapped to ChEBI in the
    # ligand_id table, use the chebiId_inchi file downloaded from the ChEBI website
    # to map them if possible.
    with zipfile.ZipFile(open(settings.BRENDA_LIGAND_ID_FNAME, 'r')) as z:
        ligand_id_table = pd.read_csv(z.open('ligand_id_table.csv'), sep=',', escapechar='\\', quotechar='"')
    
    ligand_id_table.sort_values('LigandID', inplace=True)
    
    # Add remove all references to ChEBI 33206, which is probably a mistake
    ligand_id_table.loc[ligand_id_table['chebiID'] == 'CHEBI:33206', 'chebiID'] = None
    
    # divide the table into three groups:
    # 1) ligands that are already mapped to ChEBI in the ligand_id_table
    # 2) ligands that not mapped, but have an InChI so we can try to map them
    # 3) orphan ligands, that have no ChEBI nor InChI
    mapped_ligands = ligand_id_table[~pd.isnull(ligand_id_table['chebiID'])]
    
    unmapped_ligands = ligand_id_table[pd.isnull(ligand_id_table['chebiID']) & (~pd.isnull(ligand_id_table['inchi']))]
    unmapped_ligands = unmapped_ligands.drop('chebiID', axis=1)
    
    orphan_ligands = ligand_id_table[(pd.isnull(ligand_id_table['chebiID'])) & (pd.isnull(ligand_id_table['inchi']))]
    
    print "There are %d mapped, %d unmapped and %d orphan ligands" % \
        (mapped_ligands.shape[0], unmapped_ligands.shape[0], orphan_ligands.shape[0])
        
    #%%
    
    # load the ChEBI to InChI table
    chebi_inchi_df = settings.get_chebi_inchi_df()
    chebi_inchi_df.groupby('inchi').first()
    newly_mapped_ligands = unmapped_ligands.join(chebi_inchi_df.groupby('inchi').first(), on='inchi', how='inner')
    
    print "There are %d newly mapped ligands, using direct mapping by InChI to ChEBI" % newly_mapped_ligands.shape[0]
    
    brenda2chebi = pd.concat([mapped_ligands, newly_mapped_ligands])
    brenda2chebi = brenda2chebi[['LigandID', 'chebiID']].groupby('LigandID').first()
    brenda2chebi.index.name = 'LigandID'
    
    return brenda2chebi

def rebuild_cache():
    # read the ligand ID table and remove the word "chebi:" from the beginning of the string
    #brenda2chebi = settings.get_data_df('ligand_ids')[['LigandID', 'chebiID']]
    #brenda2chebi = brenda2chebi.set_index('LigandID')
    #brenda2chebi['chebiID'] = brenda2chebi['chebiID'].apply(lambda s: s[6:] if pd.notnull(s) else pd.np.nan)
    brenda2chebi = map_brenda_to_chebi()
    
    # join the bigg IDs into the ligand table (using the ChEBIs)
    bigg = BiGG()
    kegg = KEGG()
    chebi2bigg = bigg.metabolite_df
    ligand_df = brenda2chebi.join(chebi2bigg, how='left', on='chebiID')

    # load the manually mapped ligand table
    manually_mapped_ligand_df = pd.DataFrame.from_csv(os.path.join(settings.DATA_DIR, 'ligand_ids_manual_mapping.csv'), index_col=0)
    print "There are %d manually mapped ligands" % manually_mapped_ligand_df.shape[0]
    for lid in manually_mapped_ligand_df.index:
        ligand_df.loc[lid, ['chebiID', 'bigg.metabolite']] = manually_mapped_ligand_df.loc[lid, ['chebiID', 'bigg.metabolite']]
    
    # combine the mapping from CIDs to ChEBIs with the mapping from ligand_IDs
    # to ChEBIs to get a direct mapping from BRENDA to KEGG compound IDs
    chebi2kegg = kegg.get_kegg_df()
    chebi2kegg.rename(columns={'name': 'Compound'}, inplace=True)
    chebi2kegg = chebi2kegg.groupby('chebiID').first()
    ligand_df = ligand_df.join(chebi2kegg, how='left', on='chebiID')
    
    brenda_zip = zipfile.ZipFile(open(settings.BRENDA_ZIP_FNAME, 'r'))
    for d in settings.BRENDA_INPUT:
        df = pd.DataFrame.from_csv(brenda_zip.open(d['fname'] + '.csv', 'r'), index_col=None)
        df = df.join(ligand_df, how='left', on='LigandID')
        settings.write_cache(d['fname'], df)

if __name__ == '__main__':
    rebuild_cache()