# Map activators and inhibitors to ecoli metabolic model

# 1. NB: For now, we are using all possible mappings of EC numbers, and many of them map to multiple reactions. May need to think about how to deal with this.
# 2. NB: If there is activation and inhibition, we have to deal with this carefully

import os, sys, numpy as np, pandas as pd, pdb, scipy as sp, pickle
from scipy.io import loadmat

# Load EC number dictionary
ecdict = pickle.load( open( '../cache/ecgroups.p','rb' ) )

# Read BIGG model
ecoli = loadmat('../data/iJO1366.mat')['iJO1366']
fields = ecoli.dtype.names
mnames = [ecoli['mets'][0][0][item][0][0] for item in range(len(ecoli['mets'][0][0]))]
rnames = [ecoli['rxns'][0][0][item][0][0] for item in range(len(ecoli['rxns'][0][0]))]

# Set to lower case
mnames = [unicode.lower(item) for item in mnames]
rnames = [unicode.lower(item) for item in rnames]

S = pd.DataFrame( ecoli['S'][0,0].todense(), index = mnames, columns = rnames )
R = pd.DataFrame( 0, index = mnames, columns = rnames )

# Read in the data
act = pd.read_csv('../cache/ecoli_activating_compounds_bigg.csv',header = 0,index_col = 0)
inh = pd.read_csv('../cache/ecoli_ki_bigg.csv',header = 0,index_col = 0)
vardict = {'act':act,'inh':inh}

numskipped = 0
numused = 0

for f in vardict.keys():
	
	d = vardict[ f ]
	if f == 'act':
		val = 1
	if f == 'inh':
		val = -1
	
	for row in d.index:
		
		ecnumber = d.at[ row, 'EC_number' ]
		if ecnumber in ecdict.keys():
			
			rname = ecdict[ ecnumber ]
			R.ix[ d.at[ row,'bigg.metabolite' ], rname ] += val
			numused += 1
			
		else:
			
			print('Skipping ' + ecnumber + ', not in EC dictionary.')
			numskipped += 1
			
# Write to file
R.to_csv('../cache/iJO1366_SMRN.csv')