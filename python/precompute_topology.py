# Pre-compute the shortest path length in the stoichiometric matrix

# NB: check some of the shortest path calcs?

import os, numpy as np, pdb, settings, networkx as nx, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
plt.ion()
plt.close('all')
sns.set_style('ticks')

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
metnames = [item +'_metabolite' for item in S_sign.index.tolist()]
rxnnames =  [item + '_rxn' for item in S_sign.columns.tolist()]
nodenames = metnames  + rxnnames
namedict = {item:nodenames[item] for item in range(len(nodenames))}

# Create networkx graph object
nxA = nx.to_networkx_graph( A )
nxA = nx.relabel_nodes(nxA,namedict)

# Get shortest path length and write to file
print('Calculating shortest path lengths for all pairs...')
spl = nx.shortest_path_length( nxA )

# Read in the SMRN
smrn_temp = pd.read_csv(os.path.join(settings.RESULT_DIR, 'iJO1366_SMRN.csv'),
                   header = 0,index_col = 0 )
                   
# Get unique groups
smrn = smrn_temp.groupby(['bigg.metabolite','bigg.reaction']).first().reset_index()

# Only grab data relating to edges in SMRN in order to save space
print('Grabbing only relevant shortest path lengths...')
numskipped = 0
for ii in smrn.index:

    metsource = smrn.at[ ii, 'bigg.metabolite' ]
    rxntarget = smrn.at[ ii, 'bigg.reaction' ]
    
    # We trimmed the _c, add it back
    metsource = metsource + '_c_metabolite'
    rxntarget = rxntarget+ '_rxn'
    
    if metsource not in spl.keys():
		print('Skipping this entry ' + metsource + ',' + rxntarget)
		numskipped += 1
		continue

    smrn.at[ ii, 'SPL' ] = spl[ metsource ][ rxntarget ] # divide by 2 to account for bipartite graph distance

print('Total number of entries skipped = ' + str(numskipped))

# Save data
smrn.to_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN_spl.csv'))

# Draw figure
smrn2 = smrn[pd.notnull(smrn['SPL'])]
sns.distplot( smrn2['SPL'],hist = True,kde = False, hist_kws = {'width': .5} )
plt.xlabel('Distance on Bipartite Metabolite-Enzyme Stoichiometric Matrix')
plt.ylabel('Number of SMRN Edges')
plt.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.pdf'))