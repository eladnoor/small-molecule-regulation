# Map activators and inhibitors to ecoli metabolic model

# 1. NB: For now, we are using all possible mappings of EC numbers, 
#        and many of them map to multiple reactions. 
#        May need to think about how to deal with this.
# 2. NB: If there is activation and inhibition, we have to deal with this carefully

import os
import pandas as pd
import settings

#%% Read in the data for activators and inhibitors from BRENDA
act = pd.read_csv(os.path.join(settings.CACHE_DIR, 'ecoli_activating_compounds_bigg.csv'),
                  header=0, index_col=0)
act['Value'] = 1                  
inh = pd.read_csv(os.path.join(settings.CACHE_DIR, 'ecoli_ki_bigg.csv'),
                  header=0, index_col=0)
inh['Value'] = -1
effectors = pd.concat([act, inh])
effectors.drop('Commentary', 1, inplace=True)

#%% load BiGG reaction to EC number mapping
bigg2ec = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'bigg2ec.csv'))

# Read BIGG model
model, metabolites, reactions, S = settings.get_ecoli_json()

# Set to lower case
mnames = map(unicode.lower, metabolites)
rnames = map(unicode.lower, reactions)

#%% merge BRENDA effector data with the bigg2ec table to get the BiGG reaction IDs
bigg_effectors = pd.merge(effectors, bigg2ec, how='inner', on='EC_number')

# Write to file
bigg_effectors.to_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN.csv'))

#%% group by the reaction and metabolite IDs and write to another CSV file
bigg_effectors_grouped = bigg_effectors.groupby(['bigg.reaction', 'bigg.metabolite']).sum().reset_index()
bigg_effectors_grouped.to_csv(os.path.join(settings.CACHE_DIR, 'iJO1366_SMRN_grouped.csv'))
