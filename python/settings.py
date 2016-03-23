# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:36:02 2016

@author: noore
"""

import os
import __main__ as main
import pandas as pd

SCRIPT_DIR = os.path.dirname(main.__file__)
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])
DATA_DIR = os.path.join(BASE_DIR, 'data')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')

KEGG2CHEBI_FNAME = os.path.join(CACHE_DIR, 'kegg2chebi.csv')
ECOLI_MODEL_FNAME = os.path.join(DATA_DIR, 'iJO1366.json')

def get_data_df(fname):
    return pd.DataFrame.from_csv(os.path.join(DATA_DIR, fname + '.csv'), header=0, index_col=None)
