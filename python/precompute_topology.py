# Pre-compute the shortest path length in the stoichiometric matrix

# NB: check some of the shortest path calcs?

import os, numpy as np, pdb, settings, networkx as nx, pandas as pd

# Read BIGG model
model, S = settings.get_ecoli_json()

# List of metabolites to remove from adjacency matrix
met2remove = ['h_c','h2o_c','h_p','h2o_p','co2_c','o2_c']

# Make an adjacency matrix
S_sign = S.drop( met2remove )
S_sign = np.sign( np.abs( S_sign ) )

metzero = np.zeros( (S_sign.shape[0],S_sign.shape[0]) )
rxnzero = np.zeros( (S_sign.shape[1],S_sign.shape[1]) )
A = np.bmat( [[metzero, S_sign], [S_sign.transpose(), rxnzero]] )
nodenames = S_sign.index.tolist() + S_sign.columns.tolist()
namedict = {item:nodenames[item] for item in range(len(nodenames))}

# Create networkx graph object
nxA = nx.to_networkx_graph( A )
nxA = nx.relabel_nodes(nxA,namedict)

# Get shortest path length and write to file
print('Calculating shortest path lengths for all pairs...')
spl = nx.shortest_path_length( nxA )

# Read in the SMRN
smrn = pd.read_csv( '../cache/iJO1366_SMRN.csv',header = 0,index_col = 0 )

# Drop duplicate entries if they exist
smrn = smrn.drop_duplicates()

# Only grab data relating to edges in SMRN in order to save space
print('Grabbing only relevant shortest path lengths...')
for ii in smrn.index:

    metsource = smrn.at[ ii, 'bigg.metabolite' ]
    rxntarget = smrn.at[ ii, 'bigg.reaction' ]
    smrn.at[ ii, 'SPL' ] = (spl[ metsource ][ rxntarget ] + 1) / 2 # divide by 2 to account for bipartite graph distance

# Save data
smrn.to_csv('../cache/iJO1366_SMRN_spl.csv')

