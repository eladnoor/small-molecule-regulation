# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:58:07 2016

@author: noore

This script loads the inhibition and activation links and plots them
onto the metabolic network using Cytoscape

Requirements:
sudo apt install cytoscape graphviz libgraphviz-dev libcgraph6
sudo pip install networkx pygraphviz py2cytoscape
"""

import settings as S
import pandas as pd
import os
import matplotlib.pyplot as plt
import escher

ki = pd.DataFrame.from_csv(os.path.join(S.CACHE_DIR, 'ecoli_ki_bigg.csv'))
ki_unique = ki.loc[ki.groupby(['EC_number', 'bigg.metabolite'])['Value'].idxmin()]
ki_unique.sort_values('EC_number', inplace=True)
ki_unique.to_csv('../res/ki_unique.csv')

bigg2ec = pd.DataFrame.from_csv(os.path.join(S.CACHE_DIR, 'bigg2ec.csv'))
ki_unique = pd.merge(ki_unique, bigg2ec, on='EC_number')

#act = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_activating_bigg.csv'))
#act_unique = act.groupby(['EC_number', 'bigg.metabolite']).first().reset_index()

reactions = ki_unique['bigg.reaction'].unique()
metabolites = ki_unique['bigg.metabolite'].unique()

metabolite_data = dict({m : 1 for m in metabolites})
reaction_data = dict({r.upper() : 1 for r in reactions})

my_builder = escher.Builder(map_json=os.path.join(S.DATA_DIR, 'iJO1366.CCM.json'),
                            metabolite_data=metabolite_data,
                            reaction_data=reaction_data)
my_builder.display_in_browser()

