# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:36:02 2016

@author: noore
"""

import os, sys
import __main__ as main
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(main.__file__)
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])
DATA_DIR = os.path.join(BASE_DIR, 'data')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')

KEGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'kegg2chebi.csv')
BIGG2KEGG_FNAME  = os.path.join(CACHE_DIR, 'bigg2kegg.csv')
BIGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'bigg2chebi.csv')
ECOLI_JSON_FNAME = os.path.join(DATA_DIR, 'iJO1366.json')
ECOLI_SBML_FNAME = os.path.join(DATA_DIR, 'iJO1366.xml.gz')
ECOLI_XLS_FNAME = os.path.join(DATA_DIR, 'inline-supplementary-material-2.xls')

def get_data_df(fname):
    return pd.DataFrame.from_csv(os.path.join(DATA_DIR, fname + '.csv'), header=0, index_col=None)

def get_ecoli_json():
    with open(ECOLI_JSON_FNAME) as fp:
        model = json.load(fp)
    
    # read the stoichiometric matrix from the model
    metabolites = [x['id'] for x in model['metabolites']]
    reactions = [x['id'] for x in model['reactions']]
    S = np.matrix(np.zeros((len(metabolites), len(model['reactions']))))
    for j, d in enumerate(model['reactions']):
        for met, coeff in d['metabolites'].iteritems():
            S[metabolites.index(met), j] = coeff
    
    return model, metabolites, reactions, S
    
def get_reaction_table_from_xls():
    with open(ECOLI_XLS_FNAME) as fp:
        return pd.read_excel(fp, sheetname=2, header=0)

def plotdiag():
	x1,x2,y1,y2 = plt.axis()
	minplot = np.min([x1,y1])
	maxplot = np.max([x2,y2])
	plt.plot([minplot, maxplot], [minplot,maxplot], 'k-', lw=2)
	plt.plot([minplot,maxplot/10],[minplot*10,maxplot],'k--',lw=2)
	plt.plot([minplot*10,maxplot],[minplot,maxplot/10],'k--',lw=2)
        
try:
    import cobra
    def get_ecoli_sbml():
        return cobra.io.read_sbml_model(ECOLI_SBML_FNAME)
except ImportError:
    sys.stderr.write("WARNING: please install cobrapy to have full functionality")
