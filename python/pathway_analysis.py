# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 13:55:35 2016

@author: noore

Make a bipartite graph of pathways-metabolites, with 3 types
of edges:
1) metabolite is a substrate or product of an enzyme in pathway
2) metabolite activates an enzyme in pathway
3) metabolite inhibits an enzyme in pathway
"""
from bigg import BiGG
import settings
import pandas as pd
import json
import networkx as nx
import matplotlib.pyplot as plt
from plot_figures_for_paper import FigurePlotter

bigg = BiGG()
fplot = FigurePlotter()

#%%

# make a DataFrame containing a 1-to-many mapping of 
# all Reaction IDs and their subsystems
with open(settings.ECOLI_JSON_FNAME) as fp:
    ecoli_model = json.load(fp, encoding='UTF-8')

subsystem_data = []
stoich_data = []
for r in ecoli_model['reactions']:
    rid = r['id'].lower()
    if 'subsystem' in r:
        subsystem_data.append((rid, r['subsystem']))
    if 'metabolites' in r:
        for met, coeff in r['metabolites'].iteritems():
            stoich_data.append((rid, met, coeff))

reaction_subsystem_df = pd.DataFrame(subsystem_data,
                                     columns=('bigg.reaction', 'subsystem'))
reaction_subsystem_df.set_index('bigg.reaction', inplace=True)

stoich_df = pd.DataFrame(stoich_data,
                         columns=('bigg.reaction', 'bigg.metabolite', 'coeff'))                            
          
# now associate every metabolite to subsystems by joining the two tables
                         
metabolite_subsystem_df = stoich_df.join(reaction_subsystem_df, on='bigg.reaction')
metabolite_subsystem_df.drop('bigg.reaction', axis=1, inplace=True)
metabolite_subsystem_df.drop('coeff', axis=1, inplace=True)
metabolite_subsystem_df.drop_duplicates(inplace=True)

# keep only cytoplasmic metabolites, and remove the suffix _c
metabolite_subsystem_df = metabolite_subsystem_df[metabolite_subsystem_df['bigg.metabolite'].str[-2:] == '_c']
metabolite_subsystem_df.loc[:, 'bigg.metabolite'] = metabolite_subsystem_df['bigg.metabolite'].apply(lambda s: s[0:-2].lower())

# keep only regulating metabolites
# first, we need to join the regulation table with the EC-to-bigg dataframe
reg = fplot.regulation
reg = pd.merge(reg, bigg.reaction_df, on='EC_number', how='inner')
reg = reg.join(reaction_subsystem_df, on='bigg.reaction', how='inner')

metabolite_set = set(reg['bigg.metabolite']).intersection(metabolite_subsystem_df['bigg.metabolite'])
subsystem_set = set(metabolite_subsystem_df['subsystem'].unique())
if pd.np.nan in subsystem_set:
    subsystem_set.remove(pd.np.nan)

#%%
G = nx.DiGraph()
for subsys in subsystem_set:
    G.add_node(subsys, color='blue')
for met in metabolite_set:
    G.add_node(met, color='black')
for row in metabolite_subsystem_df.iterrows():
    met = row[1]['bigg.metabolite']
    subsys = row[1]['subsystem']
    if met in metabolite_set:
        G.add_edge(met, subsys, color='black')

for row in reg.iterrows():
    met = row[1]['bigg.metabolite']
    subsys = row[1]['subsystem']
    mode = row[1]['Mode']
    if met in metabolite_set and subsys in subsystem_set:
        if mode == '+':
            G.add_edge(met, subsys, color='green')
        else:
            G.add_edge(met, subsys, color='red')

fig, ax = plt.subplots(1, 1, figsize=(20,20))
nx.draw_spring(G, ax=ax)
fig.savefig('../res/pathway_analysis.svg')