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

METS_TO_REMOVE = ['h_c', 'h2o_c', 'h_p', 'h2o_p', 'co2_c', 'o2_c']


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
    S_sparse = S_sparse[~S_sparse['bigg.metabolite'].isin(METS_TO_REMOVE)]

    B = nx.Graph()
    B.add_nodes_from(S_sparse['bigg.metabolite'].unique(), bipartite=0)
    B.add_nodes_from(S_sparse['bigg.reaction'].unique(), bipartite=1)
    B.add_weighted_edges_from(S_sparse.as_matrix())
    return B


###############################################################################
if __name__ == "__main__":
    # %% Read BIGG model
    model, S = settings.get_ecoli_json()
    B = convert_to_bipartite(S)
    spl = nx.shortest_path_length(B)
    spl_values = []
    for met in (set(S.index) - set(METS_TO_REMOVE)):
        rxns = set(S.columns).intersection(spl[met].keys())
        spl_values += map(spl[met].get, rxns)
    spl_values = (np.array(spl_values) - 1.0) / 2

    smrn = pd.read_csv(os.path.join(settings.CACHE_DIR,
                                    'iJO1366_SMRN.csv'), index_col=None)
    smrn_dist = smrn[['bigg.metabolite', 'bigg.reaction']].drop_duplicates()
    smrn_dist['distance'] = pd.np.nan
    for i, row in smrn_dist.iterrows():
        source = row['bigg.metabolite'] + '_c'  # remember we dropped it before
        target = row['bigg.reaction']
        if source in METS_TO_REMOVE:
            continue
        smrn_dist.at[i, 'distance'] = (spl[source][target] - 1.0) / 2.0

    # %% Save data
    smrn_dist = smrn_dist.dropna()
    smrn_dist.to_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN_dist.csv'),
                     index=False)

    # %% Draw figure
    Nmax = 7
    fig, axs = plt.subplots(2, 1, figsize=(6, 8), sharex=False)
    axs[0].hist(smrn_dist['distance'], alpha=1, normed=True, align='left',
                bins=np.arange(Nmax+1), rwidth=.8, label='SMRN pairs',
                color='#8a5ec9', linewidth=0)
    axs[1].hist(spl_values, alpha=1, normed=True, align='left',
                bins=np.arange(Nmax+1), rwidth=.8, label='All pairs',
                color='#939598', linewidth=0)

    axs[0].set_xticks(np.arange(Nmax))
    axs[1].set_xticks(np.arange(Nmax))
    axs[0].set_xlim(-1, Nmax)
    axs[1].set_xlim(-1, Nmax)
    axs[0].set_ylabel('Fraction of metabolite-enzyme pairs')
    axs[0].set_xlabel('')
    axs[1].set_ylabel('Fraction of metabolite-enzyme pairs')
    axs[1].set_xlabel('Distance in # reactions between metabolite and enzyme')

    axs[0].annotate('a', xy=(0.02, 0.98),
                    xycoords='axes fraction', ha='left', va='top',
                    size=20)
    axs[0].annotate('SMRN pairs', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    axs[1].annotate('b', xy=(0.02, 0.98),
                    xycoords='axes fraction', ha='left', va='top',
                    size=20)
    axs[1].annotate('All pairs', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    plt.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.pdf'))

    # %% Draw figure where data is split between modes
    smrn_merged = pd.merge(smrn, smrn_dist, on=['bigg.metabolite',
                                                'bigg.reaction'])
    dist_mode_df = smrn_merged.groupby(('bigg.metabolite',
                                        'bigg.reaction', 'Mode')).first()
    dist_mode_df = dist_mode_df[['distance']].reset_index()
    inh_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '-', 'distance']
    act_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '+', 'distance']

    # Make a histogram indicating the relative distance
    # activating/inhibiting edges travel
    Nmax = 7
    fig, axs = plt.subplots(2, 1, figsize=(6, 8), sharex=False)
    axs[0].hist(inh_distances, alpha=1, normed=True, align='left',
                bins=np.arange(Nmax+1), rwidth=.8, label='Inhibition',
                color='#f5821f', linewidth=0)
    axs[1].hist(act_distances, alpha=1, normed=True, align='left',
                bins=np.arange(Nmax+1), rwidth=.8, label='Activation',
                color='#00aeef', linewidth=0)

    axs[0].set_xticks(np.arange(Nmax))
    axs[1].set_xticks(np.arange(Nmax))
    axs[0].set_xlim(-1, Nmax)
    axs[1].set_xlim(-1, Nmax)
    axs[0].set_ylabel('Fraction of metabolite-enzyme pairs')
    axs[0].set_xlabel('')
    axs[1].set_ylabel('Fraction of metabolite-enzyme pairs')
    axs[1].set_xlabel('Distance in # reactions between metabolite and enzyme')

    axs[0].annotate('a', xy=(0.02, 0.98),
                    xycoords='axes fraction', ha='left', va='top',
                    size=20)
    axs[0].annotate('Inhibition', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    axs[1].annotate('b', xy=(0.02, 0.98),
                    xycoords='axes fraction', ha='left', va='top',
                    size=20)
    axs[1].annotate('Activation', xy=(0.9, 0.9),
                    xycoords='axes fraction', ha='right', va='top',
                    size=14)
    plt.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_mode_distances.pdf'))
