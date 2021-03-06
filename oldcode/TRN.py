# Script to compare TRN and SMRN in E. Coli

import os, sys, settings, pandas as pd, matplotlib.pyplot as plt, pdb, numpy as np, seaborn as sns
plt.ion()
plt.close('all')
sns.set_style('ticks')

# Read BIGG model
model, S = settings.get_ecoli_json()

# Get the unique genes and reactions
genes = model['genes']
rxns = model['reactions']
rxnnames = [item['id'] for item in rxns]

# Dictionary mapping b-numbers to gene names
gdict = {item['id']:item['name'] for item in genes}

# Make a dataframe indicating regulatory relationships
r2g = pd.DataFrame( 0,columns = gdict.values(), index = rxnnames )

# Parse rxns to get all genes regulating that reaction
for rxn in rxns:
    
    rxnname = rxn['id']
    
    if rxn['gene_reaction_rule'] == '':
        continue
        
    bsplit = rxn['gene_reaction_rule'].split('b')
    bIDs = ['b' + item[:4] for item in bsplit if item != '']
    gIDs = [gdict[item] for item in bIDs if item in gdict.keys()]
    
    r2g.ix[ rxnname, gIDs ] = 1

# Make index of r2g lowercase
r2g.index = [item.lower() for item in r2g.index]

# Get TRN
trn = pd.read_csv('../data/RegulonDB_network_tf_gene.txt',header = None,comment = '#',sep = '\t')
trn.columns = ['tf','gene','effect','evidence','level','none']

# Drop any entries where the strength of the evidence is NA (i.e. not "Weak" or "Strong")
trn = trn[ pd.notnull(trn['level']) ]

# For each gene, count the number of entries
tfcounts = trn.groupby('gene').size()


# Get number of genes regulating each reaction, and number of TFs
gsums = r2g.sum(axis = 1)

ixgenes = np.intersect1d( tfcounts.index, r2g.columns )
tfsums = r2g.ix[:,ixgenes].dot(tfcounts.ix[ixgenes])

# Read in the SMRN
smrn = pd.read_csv( '../cache/iJO1366_SMRN.csv',header = 0,index_col = 0 )

# Drop duplicate entries if they exist
smrn = smrn.drop_duplicates()

# Get some statistics on reactions in the SMRN
smrn_rxnsums = smrn.groupby('bigg.reaction').size()
smrn_actsums = smrn[ smrn['Value'] == 1 ].groupby('bigg.reaction').size()
smrn_inhsums = smrn[ smrn['Value'] == -1].groupby('bigg.reaction').size()

alld = pd.concat( [smrn_rxnsums, gsums, tfsums], join = 'inner',axis = 1 )
alld.columns = ['SMRNAll','Genes','TFs']

# Add information on just activating/inhibiting edges in SMRN. These are not included in pd.concat above because they might lack some indices.
alld.ix[ smrn_actsums.index, 'SMRNAct' ] = smrn_actsums
alld.ix[ smrn_inhsums.index, 'SMRNInh' ] = smrn_inhsums

# Fill in spots that had no activating/inhibiting entries
alld = alld.fillna(0)

# Make some pivot tables
smrn_vs_genes = pd.pivot_table(alld,index=['SMRNAll', 'Genes'],aggfunc=len).reset_index()

smrn_vs_tfs = pd.pivot_table(alld, index=['SMRNAll', 'TFs'], aggfunc=len).reset_index()

# Plot result
plt.figure(figsize=(10, 10))
plt.subplot(2,1,1)
plt.scatter (smrn_vs_genes['SMRNAll'], smrn_vs_genes['Genes'], s=5*smrn_vs_genes['SMRNAct'] )
plt.xlabel('# of Regulators for Reaction')
plt.ylabel('# of Genes Encoding Reaction')

plt.subplot(2,1,2)
plt.scatter (smrn_vs_tfs['SMRNAll'], smrn_vs_tfs['TFs'], s=5*smrn_vs_tfs['SMRNAct'] )
plt.xlabel('# of Regulators for Reaction')
plt.ylabel('# of TFs Regulating Reaction')

plt.savefig( '../res/TRN_vs_SMRN.svg' )