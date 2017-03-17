# Pre-compute the shortest path length in the stoichiometric matrix

# NB: check some of the shortest path calcs?

import pdb
import settings
import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.ion()
plt.close('all')
sns.set_style('ticks')

METS_TO_REMOVE = ['h', 'h2o', 'co2', 'o2', 'pi', 'atp', 'adp', 'amp',
                  'nad', 'nadh', 'nadp', 'nadph', 'coa', 'thf', '5mthf',
                  '5fthf', 'methf', 'mlthf', 'nh4', 'cmp', 'q8', 'q8h2',
                  'udp', 'udpg', 'fad', 'fadh2', 'ade', 'ctp', 'gtp', 'h2o2',
                  'mql8', 'mqn8', 'na1', 'ppi', 'acp']


def convert_to_bipartite(S):
    """
        convert a standard stoichiometric matrix (in a Pandas DataFrame)
        to a bipartite graph with an edge between every reactant and all its
        reactions
    """
    # convert the stoichiometric matrix to a sparse representation
    S_sparse = pd.melt(S.reset_index(),
                       id_vars='bigg.metabolite', value_name='coeff')
    S_sparse = S_sparse[S_sparse.coeff != 0]

    # remove the high-degree metabolites that we want to ignore for graph
    # distance
    met_comp = S_sparse['bigg.metabolite'].str.rsplit('_', 1, expand=True)
    S_sparse = S_sparse[(~met_comp[0].isin(METS_TO_REMOVE)) & (met_comp[1] == 'c')]
    S_sparse['bigg.metabolite'] = met_comp[0].str.upper()

    mets = set(S_sparse['bigg.metabolite'].unique())
    rxns = set(S_sparse['bigg.reaction'].unique())

    B = nx.Graph()
    B.add_nodes_from(mets, bipartite=0)
    B.add_nodes_from(rxns, bipartite=1)
    B.add_weighted_edges_from(S_sparse.as_matrix())
    return B, mets, rxns


def calculate_distances(smrn):

    smrn['bigg.metabolite'] = smrn['bigg.metabolite'].str.upper()

    # %% Read BIGG model
    model, S = settings.get_ecoli_json()
    B, mets, rxns = convert_to_bipartite(S)
    spl = nx.shortest_path_length(B)
    spl_values = []
    for met in mets:
        spl_values += map(spl[met].get, rxns.intersection(spl[met].keys()))
    all_distances = (np.array(spl_values) - 1.0) / 2

    smrn_dist = smrn[['bigg.metabolite', 'bigg.reaction']].drop_duplicates()
    smrn_dist['distance'] = pd.np.nan
    for i, row in smrn_dist.iterrows():
        source = row['bigg.metabolite']  # remember we dropped it before
        target = row['bigg.reaction']
        if source.lower() in METS_TO_REMOVE:
            continue
        if target in spl[source]:
            smrn_dist.at[i, 'distance'] = (spl[source][target] - 1.0) / 2.0

    # %% Save data
    smrn_dist = smrn_dist.dropna()
    return smrn_dist, all_distances
