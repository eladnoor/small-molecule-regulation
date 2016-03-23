# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 23:42:54 2016

@author: eladn
"""

import json
import settings
import numpy as np

with open(settings.ECOLI_MODEL_FNAME) as fp:
    model = json.load(fp)

# read the stoichiometric matrix from the model
metabolites = [x['id'] for x in model['metabolites']]
S = np.matrix(np.zeros((len(metabolites), len(model['reactions']))))
for j, d in enumerate(model['reactions']):
    for met, coeff in d['metabolites'].iteritems():
        S[metabolites.index(met), j] = coeff

