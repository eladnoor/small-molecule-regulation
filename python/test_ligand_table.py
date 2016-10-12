# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 18:40:38 2016

@author: noore
"""
import settings
import pandas as pd
import zipfile

# For all the Ligands that have an InChI but are not mapped to ChEBI in the
# ligand_id table, use the chebiId_inchi file downloaded from the ChEBI website
# to map them if possible.
with zipfile.ZipFile(open(settings.BRENDA_LIGAND_ID_FNAME, 'r')) as z:
    ligand_id_table = pd.read_csv(z.open('ligand_id_table.csv'), sep=',', escapechar='\\', quotechar='"')

ligand_id_table.sort_values('LigandID', inplace=True)

# Add NAD(P)+ to the mapping manually
ligand_id_table.loc[ligand_id_table['LigandID'] == 7, 'chebiID'] = 'CHEBI:15846'  # NAD+
ligand_id_table.loc[ligand_id_table['LigandID'] == 10, 'chebiID'] = 'CHEBI:18009' # NADP+
ligand_id_table.loc[ligand_id_table['LigandID'] == 161, 'chebiID'] = 'CHEBI:33191' # KCN
ligand_id_table.loc[ligand_id_table['chebiID'] == 'CHEBI:33206', 'chebiID'] = None

# divide the table into three groups:
# 1) ligands that are already mapped to ChEBI in the ligand_id_table
# 2) ligands that not mapped, but have an InChI so we can try to map them
# 3) orphan ligands, that have no ChEBI nor InChI
mapped_ligands = ligand_id_table[~pd.isnull(ligand_id_table['chebiID'])]
mapped_ligands = mapped_ligands[['LigandID', 'chebiID']].groupby('LigandID').first()

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
newly_mapped_ligands = newly_mapped_ligands[['LigandID', 'chebiID']].groupby('LigandID').first()

print "There are %d newly mapped ligands, using direct mapping by InChI to ChEBI" % newly_mapped_ligands.shape[0]

all_mapped = pd.concat([mapped_ligands, newly_mapped_ligands])
all_mapped.index.name = 'LigandID'
