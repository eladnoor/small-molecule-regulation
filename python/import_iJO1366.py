# Import iJO1366 from .mat and format for mapping of regulatory interactions from BRENDA

import os
import pandas as pd
import pdb
import settings

# Read BIGG model
model, metabolites, reactions, S = settings.get_ecoli_json()

# Set to lower case
rnames = map(unicode.lower, reactions)

# Read mapping from KEGG IDs
model_reactions = pd.read_excel('../data/inline-supplementary-material-2.xls', sheetname=2, header=0)
bigg2ec = model_reactions[['Reaction Abbreviation', 'EC Number']]
bigg2ec.rename(columns={'Reaction Abbreviation': 'bigg.reaction',
                        'EC Number': 'EC_number'},
                        inplace=True)
bigg2ec = bigg2ec.loc[~bigg2ec['EC_number'].isnull()]

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
	print('Multiple reactions with the same name!')
	pdb.set_trace()

# write bigg2ec mapping to csv file
bigg2ec.to_csv(os.path.join(settings.CACHE_DIR, 'bigg2ec.csv'))
