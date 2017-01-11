# Pre-compute the shortest path length in the stoichiometric matrix

# NB: check some of the shortest path calcs?

import os
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

    # remove the high-degree metabolites that are we want to ignore for graph
    # distance
    met_suff = S_sparse['bigg.metabolite'].str.rsplit('_', 1, expand=True)
    S_sparse['bigg.metabolite'] = met_suff[0]
    S_sparse['compartment'] = met_suff[1]
    S_sparse = S_sparse[~S_sparse['bigg.metabolite'].isin(METS_TO_REMOVE)]
    S_sparse = S_sparse[S_sparse['compartment'] == 'c']
    S_sparse.drop('compartment', axis=1, inplace=True)

    mets = set(S_sparse['bigg.metabolite'].unique())
    rxns = set(S_sparse['bigg.reaction'].unique())

    B = nx.Graph()
    B.add_nodes_from(mets, bipartite=0)
    B.add_nodes_from(rxns, bipartite=1)
    B.add_weighted_edges_from(S_sparse.as_matrix())
    return B, mets, rxns


###############################################################################
if __name__ == "__main__":
    # %% Read BIGG model
    model, S = settings.get_ecoli_json()
    B, mets, rxns = convert_to_bipartite(S)
    spl = nx.shortest_path_length(B)
    spl_values = []
    for met in mets:
        spl_values += map(spl[met].get, rxns.intersection(spl[met].keys()))
    spl_values = (np.array(spl_values) - 1.0) / 2

    smrn = pd.read_csv(os.path.join(settings.CACHE_DIR,
                                    'iJO1366_SMRN.csv'), index_col=None)
    smrn_dist = smrn[['bigg.metabolite', 'bigg.reaction']].drop_duplicates()
    smrn_dist['distance'] = pd.np.nan
    for i, row in smrn_dist.iterrows():
        source = row['bigg.metabolite']  # remember we dropped it before
        target = row['bigg.reaction']
        if source in METS_TO_REMOVE:
            continue
        if target in spl[source]:
            smrn_dist.at[i, 'distance'] = (spl[source][target] - 1.0) / 2.0

    # %% Save data
    smrn_dist = smrn_dist.dropna()
    smrn_dist.to_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN_dist.csv'),
                     index=False)

    # %% Draw figure
    smrn_merged = pd.merge(smrn, smrn_dist, on=['bigg.metabolite',
                                                'bigg.reaction'])
    dist_mode_df = smrn_merged.groupby(('bigg.metabolite',
                                        'bigg.reaction', 'Mode')).first()
    dist_mode_df = dist_mode_df[['distance']].reset_index()
    inh_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '-', 'distance']
    act_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '+', 'distance']

    Nmax = 9
    bins = range(Nmax+1)
    args = {'alpha': 1, 'normed': True, 'align': 'left', 'bins': bins,
            'linewidth': 0, 'rwidth': 0.8}
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=False)
    axs[0, 0].hist(spl_values, color='#939598', **args)
    axs[0, 1].hist(smrn_dist['distance'], color='#8a5ec9', **args)
    axs[1, 0].hist(inh_distances, color='#f5821f', **args)
    axs[1, 1].hist(act_distances, color='#00aeef', **args)

    for i in range(2):
        axs[i, 0].set_ylabel('Fraction of metabolite-enzyme pairs')
        axs[1, i].set_xlabel('Distance in # reactions between metabolite and enzyme')
    for i, ax in enumerate(axs.flat):
        ax.annotate(chr(ord('a') + i), xy=(0.02, 0.98),
                    xycoords='axes fraction', ha='left', va='top',
                    size=20)
        ax.set_xticks(np.arange(Nmax+1))
        ax.set_xlim(-1, Nmax+1)

    axs[0, 0].annotate('all iJO1366 pairs', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    axs[0, 1].annotate('all SMRN pairs', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    axs[1, 0].annotate('only inhibition', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    axs[1, 1].annotate('only activation', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.pdf'))
    fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.png'))
