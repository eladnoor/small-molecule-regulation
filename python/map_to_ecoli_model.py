# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 23:42:54 2016

@author: eladn
"""

import settings
import os
import pandas as pd

model, metabolites, reactions, S = settings.get_ecoli_json()
ki = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_ki_bigg.csv'))
activators = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_activating_compounds_bigg.csv'))
