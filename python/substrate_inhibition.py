#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 13:11:45 2017

@author: noore
"""
import settings
import os
import pandas as pd
import json

# convert S into a table
with open(settings.ECOLI_JSON_FNAME) as fp:
    model = json.load(fp)

sparse = []
for reaction in model['reactions']:
    for met, coeff in reaction['metabolites'].iteritems():
        sparse.append([reaction['id'].lower(), met.lower(), coeff])

S_sparse = pd.DataFrame(sparse, columns=['bigg.reaction',
                                         'bigg.metabolite', 'stoichiometry'])

# remove the high-degree metabolites that we want to ignore for graph
# distance
met_comp = S_sparse['bigg.metabolite'].str.rsplit('_', 1, expand=True)
S_sparse['bigg.metabolite'] = met_comp[0].str.upper()

smrn = pd.read_csv(os.path.join(settings.CACHE_DIR,
                                'iJO1366_SMRN.csv'), index_col=None)

smrn['bigg.metabolite'] = smrn['bigg.metabolite'].str.upper()

smrn = smrn[smrn['Mode'] == '-']
subs_inhibition = pd.merge(smrn, S_sparse, on=['bigg.metabolite', 'bigg.reaction'])
subs_inhibition = subs_inhibition[subs_inhibition['stoichiometry'] < 0]

subs_inhibition.to_csv(os.path.join(settings.RESULT_DIR, 'reactant_inhibition_full.csv'))

subs_inhibition = subs_inhibition[['EC_number', 'bigg.metabolite', 'bigg.reaction']].drop_duplicates()
subs_inhibition.to_excel(os.path.join(settings.RESULT_DIR, 'reactant_inhibition.xls'))
