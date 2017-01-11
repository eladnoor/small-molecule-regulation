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
    # ligand_id table, use the chebiId_inchi file downloaded from the ChEBI
    # website to map them if possible.
    with zipfile.ZipFile(open(settings.BRENDA_LIGAND_ID_FNAME, 'r')) as z:
        ligand_id_table = pd.read_csv(z.open('ligand_id_table.csv'), sep=',',
                                      escapechar='\\', quotechar='"')

    ligand_id_names = ligand_id_table[['LigandID', 'LigandName']]
    ligand_id_table = ligand_id_table.groupby('LigandID').first().reset_index()

    ligand_id_table.sort_values('LigandID', inplace=True)

    # Add remove all references to CHEBI:33206, which is probably a mistake
    ligand_id_table.loc[ligand_id_table['chebiID'] == 'CHEBI:33206', 'chebiID'] = \
        None

    # map beta-fructose-1,6P instances to the more canonical ID used by BiGG
    ligand_id_table.loc[ligand_id_table['chebiID'] == 'CHEBI:28013', 'chebiID'] = \
        'CHEBI:37736'

    # replace the ID for glycerol-1-phosphate with glycerol-3-phosphate which
    # is much more common in models (CHEBI:14336 -> CHEBI:57597)
    ligand_id_table.loc[ligand_id_table['chebiID'] == 'CHEBI:14336', 'chebiID'] = \
        'CHEBI:57597'

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

    # load the ChEBI to InChI table
    chebi_inchi_df = settings.get_chebi_inchi_df()
    chebi_inchi_df.groupby('inchi').first()
    newly_mapped_ligands = unmapped_ligands.join(chebi_inchi_df.groupby('inchi').first(), on='inchi', how='inner')

    print "There are %d newly mapped ligands, using direct mapping by InChI to ChEBI" % newly_mapped_ligands.shape[0]

    brenda2chebi = pd.concat([mapped_ligands, newly_mapped_ligands])
    brenda2chebi = brenda2chebi[['LigandID', 'chebiID']].groupby('LigandID').first()
    brenda2chebi = ligand_id_names.join(brenda2chebi, on='LigandID', how='left')
    brenda2chebi.set_index('LigandID', inplace=True)
    brenda2chebi.index.name = 'LigandID'

    # join the bigg IDs into the ligand table (using the ChEBIs)
    bigg = BiGG()
    kegg = KEGG()
    chebi2bigg = bigg.metabolite_df
    ligand_df = brenda2chebi.join(chebi2bigg, on='chebiID', how='left')

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

    return chebi2kegg, chebi2bigg, ligand_df


def rebuild_cache():

    # read the ligand ID table and remove the word "chebi:"
    # from the beginning of the string
    chebi2kegg, chebi2bigg, ligand_df = map_brenda_to_chebi()

    # map the BRENDA data that has nothing to do with
    # regulation (i.e. MM-kinetics)
    brenda_zip = zipfile.ZipFile(open(settings.BRENDA_ZIP_FNAME, 'r'))
    for d in ['km', 'turnover']:
        df = pd.DataFrame.from_csv(brenda_zip.open(d + '.csv', 'r'),
                                   index_col=None)
        df = df.join(ligand_df, how='left', on='LigandID')
        settings.write_cache(d, df)

    act_df = pd.DataFrame.from_csv(brenda_zip.open('activating.csv', 'r'),
                                   index_col=None)
    inh_df = pd.DataFrame.from_csv(brenda_zip.open('inhibiting.csv', 'r'),
                                   index_col=None)
    ki_df = pd.DataFrame.from_csv(brenda_zip.open('ki.csv', 'r'),
                                  index_col=None)
    act_df['Mode'] = '+'
    ki_df['Mode'] = '-'
    inh_df['Mode'] = '-'
    reg_df = pd.concat([act_df, inh_df, ki_df], ignore_index=True)

    # we also need to map the ecocyc data using chebis
    reg_df = reg_df.join(ligand_df, how='left', on='LigandID')
    reg_df['Source'] = 'BRENDA'
    reg_df['Commentary'].fillna('', inplace=True)

    # load the EcoCyc data and merge with BRENDA
    ecocyc_reg_df = pd.DataFrame.from_csv(settings.ECOCYC_REG_FNAME,
                                          index_col=None)
    # replace the EcoCyc chebi ID for fructose-1,6P to be more consistent with
    # BiGG (it is the same compound, except in EcoCyc they use the -4 charges
    # pseudoisomer)
    ecocyc_reg_df.loc[ecocyc_reg_df['chebiID'] == 'CHEBI:32966', 'chebiID'] = \
        'CHEBI:37736'
    ecocyc_reg_df['Commentary'] = ecocyc_reg_df['ReactionID'].map(
        lambda i: 'ReactionID=%s' % i)
    ecocyc_reg_df['Source'] = 'EcoCyc'
    # convert uM to mM
    ecocyc_reg_df['KI_Value'] = ecocyc_reg_df['KI_Value_uM'] * 1e-3
    ecocyc_reg_df.drop('KI_Value_uM', axis=1, inplace=True)
    ecocyc_reg_df.drop('ReactionID', axis=1, inplace=True)
    ecocyc_reg_df = ecocyc_reg_df.join(chebi2bigg, how='left', on='chebiID')
    ecocyc_reg_df = ecocyc_reg_df.join(chebi2kegg, how='left', on='chebiID')

    # load the manually curated KI data
    manual_df = pd.DataFrame.from_csv(settings.MANUAL_REG_FNAME,
                                      index_col=None)
    reg_df = pd.concat([reg_df, ecocyc_reg_df, manual_df], ignore_index=True)

    settings.write_cache('regulation', reg_df)
    return reg_df

if __name__ == '__main__':
    reg_df = rebuild_cache()

    k = reg_df.copy()
    ind_with_ki = k[~pd.isnull(k['KI_Value'])].index
    k.loc[ind_with_ki, 'Mode'] = 'ki'
    k = k[['Mode', 'bigg.metabolite', 'EC_number', 'Commentary', 'Organism']]
    k['Counter'] = 1

    print "All Entries"
    print k.groupby('Mode').sum()['Counter']
    k = k[k['Organism'] == 'Escherichia coli']

    print "Only E. coli"
    print k.groupby('Mode').sum()['Counter']

    # filter out mutated enzymes
    k = k[(pd.isnull(k['Commentary'])) |
          ((k['Commentary'].str.find('mutant') == -1) &
           (k['Commentary'].str.find('mutation') == -1))]

    print "filtering mutants"
    print k.groupby('Mode').sum()['Counter']

    # remove values with unmatched ligand
    bigg = BiGG()
    k = k[pd.notnull(k['bigg.metabolite'])]
    k['bigg.metabolite'] = k['bigg.metabolite'].str.lower()
    k = k[k['EC_number'].isin(bigg.get_native_EC_numbers())]
    print "leaving only EC numbers from BiGG"
    print k.groupby('Mode').sum()['Counter']

    k = k.groupby(('bigg.metabolite', 'EC_number', 'Mode')).first().reset_index()
    print "Unique met-EC pairs"
    print k.groupby('Mode').sum()['Counter']

    print "Unique metabolites"
    print k.groupby(('bigg.metabolite', 'Mode')).first().reset_index().groupby('Mode').sum()['Counter']

    print "Unique EC numbers"
    print k.groupby(('EC_number', 'Mode')).first().reset_index().groupby('Mode').sum()['Counter']
