# Import iJO1366 from .mat and format for mapping of regulatory interactions from BRENDA

import os, sys, numpy as np, pandas as pd, pdb, scipy as sp, pickle
from scipy.io import loadmat

# Read BIGG model
ecoli = loadmat('../data/iJO1366.mat')['iJO1366']
fields = ecoli.dtype.names
rnames = [ecoli['rxns'][0][0][item][0][0] for item in range(len(ecoli['rxns'][0][0]))]

# Set to lower case
rnames = [unicode.lower(item) for item in rnames]

# Read mapping from KEGG IDs
model_reactions = pd.read_excel('../data/inline-supplementary-material-2.xls', sheetname=2, header=0)
bigg2ec = model_reactions[['Reaction Abbreviation', 'EC Number']]
bigg2ec.rename(columns={'Reaction Abbreviation': 'bigg.reaction'}, inplace=True)
bigg2ec = bigg2ec.loc[~bigg2ec['EC Number'].isnull()]

# Change all reaction IDs to lower-case (apparently the standards have changed
# since the model was published, and cases are different now).

bigg2ec['bigg.reaction'] = bigg2ec['bigg.reaction'].apply(unicode.lower)

# Confirm that the name of the reaction is the same in the .xls and the .mat file
error_rnames = [item for item in bigg2ec.index if bigg2ec.at[ item, 'bigg.reaction' ] != rnames[item] ]

if len(error_rnames) != 0:
	print('These reactions may be mapped incorrectly, proceed with care. We will replace the names in bigg2ec for consistency:')
	
	xls_values = bigg2ec.loc[ error_rnames, 'bigg.reaction' ] # Note that we use .loc to indicate that this is the proper name of the index, not its numerical order
	mat_values = [rnames[item] for item in error_rnames]
	
	print xls_values
	print mat_values
	
	bigg2ec.loc[ error_rnames, 'bigg.reaction' ] = mat_values
	
# Confirm that there are no duplicate bigg.reaction names, and if so, use this as the index
if bigg2ec['bigg.reaction'].duplicated().any():
	print( 'Multiple reactions with the same name!')
	pdb.set_trace()
else:
	bigg2ec.index = bigg2ec['bigg.reaction']

# Make a dictionary which maps from EC number to unique reaction name
ecgroups = bigg2ec.groupby('EC Number')


# Save the result
pickle.dump( ecgroups.groups, open( '../cache/ecgroups.p', 'wb' ) )