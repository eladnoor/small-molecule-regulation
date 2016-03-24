# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:13:37 2016

@author: noore
"""

import settings
import pandas as pd
import os

model = settings.get_ecoli_sbml()

bigg2kegg = []
for m in model.metabolites:
    kegg_id_list = m.annotation.get('kegg.compound')
    if kegg_id_list is None:
        kegg_id_list = []
    elif type(kegg_id_list) == str:
        kegg_id_list = [kegg_id_list]
    
    for kegg_id in kegg_id_list:
        if kegg_id[0] == 'C':
            bigg2kegg.append((m.id, kegg_id))

bigg2kegg = pd.DataFrame(bigg2kegg, columns=['bigg.metabolite', 'kegg.compound'])
bigg2kegg.to_csv(os.path.join(settings.CACHE_DIR, settings.BIGG2KEGG_FNAME))

bigg2chebi = []
for m in model.metabolites:
    chebi_id_list = m.annotation.get('chebi')
    if chebi_id_list is None:
        chebi_id_list = []
    elif type(chebi_id_list) == str:
        chebi_id_list = [chebi_id_list]
    
    for chebi_id in chebi_id_list:
        bigg2chebi.append((m.id, chebi_id))

bigg2chebi = pd.DataFrame(bigg2chebi, columns=['bigg.metabolite', 'chebi'])
bigg2chebi.to_csv(os.path.join(settings.CACHE_DIR, settings.BIGG2CHEBI_FNAME))

