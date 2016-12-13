# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:36:02 2016

@author: noore
"""

import os, sys, re
import __main__ as main
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import imp
import json
import inspect
import urllib2
from contextlib import closing

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#SCRIPT_DIR = os.path.dirname(os.path.abspath(main.__file__))
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])

DATA_DIR = os.path.join(BASE_DIR, 'data')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')
RESULT_DIR = os.path.join(BASE_DIR, 'res')

CHEBI2INCHI_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv'
KEGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'kegg2chebi.csv')
BIGG_METABOLITE_FNAME = os.path.join(DATA_DIR, 'bigg_models_metabolites.txt')
BIGG_REACTION_FNAME = os.path.join(DATA_DIR, 'bigg_models_reactions.txt')
BIGG2KEGG_FNAME  = os.path.join(CACHE_DIR, 'bigg2kegg.csv')
BIGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'bigg2chebi.csv')
BIGG2EC_FNAME = os.path.join(CACHE_DIR, 'bigg2ec.csv')
ECOLI_THERMO_CACHE_FNAME = os.path.join(CACHE_DIR, 'ecoli_thermodynamics.csv')

BRENDA_ZIP_FNAME = os.path.join(DATA_DIR, 'brenda_query_2016_06_20.zip')
BRENDA_LIGAND_ID_FNAME = os.path.join(DATA_DIR, 'ligand_id_table.csv.zip')
ECOLI_JSON_FNAME = os.path.join(DATA_DIR, 'iJO1366.json')
ECOLI_SBML_FNAME = os.path.join(DATA_DIR, 'iJO1366.xml')
ECOLI_MODEL_METABOLITES = os.path.join(DATA_DIR, 'iJO1366_model_metabolites.csv')
ECOLI_MODEL_REACTIONS = os.path.join(DATA_DIR, 'iJO1366_model_reactions.csv')
ECOLI_XLS_FNAME = os.path.join(DATA_DIR, 'inline-supplementary-material-2.xls')
ECOLI_METAB_FNAME = os.path.join(DATA_DIR, 'ecoli_metabolites_gerosa2015.csv')
ECOLI_PROT_FNAME = os.path.join(DATA_DIR, 'ecoli_proteins_schmidt2015.csv')
ECOLI_FLUX_FNAME = os.path.join(DATA_DIR, 'ecoli_fluxes_fba.csv')
ECOLI_CCM_FNAME = os.path.join(DATA_DIR, 'ecoli_ccm_ec_numbers.csv')

BRENDA_INPUT = [{'fname': 'turnover',   'value_col': 'Turnover_Number'},
                {'fname': 'ki',         'value_col': 'KI_Value'},
                {'fname': 'km',         'value_col': 'KM_Value'},
                {'fname': 'activating', 'value_col': None},
                {'fname': 'inhibiting', 'value_col': None},
                {'fname': 'CCM_Reactions', 'value_col': None} ]

def get_data_df(fname):
    return pd.DataFrame.from_csv(os.path.join(DATA_DIR, fname + '.csv'), header=0, index_col=None)

def write_cache(fname, df):
    df.to_csv(os.path.join(CACHE_DIR, fname + '.csv'))

def read_cache(fname):
    return pd.DataFrame.from_csv(os.path.join(CACHE_DIR, fname + '.csv'))

def get_reaction_table_from_xls():
    with open(ECOLI_XLS_FNAME) as fp:
        return pd.read_excel(fp, sheetname=2, header=0)

def plotdiag(lw=2, ax=None):
    if ax is None:
        ax = plt
    x1, x2, y1, y2 = ax.axis()
    minplot = np.min([x1, y1])
    maxplot = np.max([x2, y2])
    ax.plot([minplot,    maxplot],    [minplot,    maxplot],    'k-',  lw=lw)
    ax.plot([minplot,    maxplot/10], [minplot*10, maxplot],    'k--', lw=lw)
    ax.plot([minplot*10, maxplot],    [minplot,    maxplot/10], 'k--', lw=lw)
        
try:
    imp.find_module('cobra')
    cobrafound = True
except ImportError:
    cobrafound = False
    sys.stderr.write("WARNING: please install cobrapy to have full functionality")
if cobrafound:
    def get_ecoli_sbml():
        return cobra.io.read_sbml_model(ECOLI_SBML_FNAME)    

def get_ecoli_json():
    with open(ECOLI_JSON_FNAME) as fp:
        model = json.load(fp)

    sparse = []
    for reaction in model['reactions']:
        for met, coeff in reaction['metabolites'].iteritems():
            sparse.append([reaction['id'].lower(), met.lower(), coeff])
    
    sparse = pd.DataFrame(sparse, columns=['bigg.reaction', 'bigg.metabolite', 'stoichiometry'])
    S = sparse.pivot(index='bigg.metabolite', columns='bigg.reaction', values='stoichiometry')
    S.fillna(0, inplace=True)
    return model, S
    
def get_chebi_inchi_df():
    with closing(urllib2.urlopen(CHEBI2INCHI_URL, 'file')) as r:
        chebi_inchi_df = pd.read_csv(r, sep='\t')
    chebi_inchi_df['chebiID'] = chebi_inchi_df['CHEBI_ID'].apply(lambda c: 'CHEBI:%d' % c)
    chebi_inchi_df.rename(columns={'InChI':'inchi'}, inplace=True)
    chebi_inchi_df.set_index('CHEBI_ID', inplace=True)
    return chebi_inchi_df
    
