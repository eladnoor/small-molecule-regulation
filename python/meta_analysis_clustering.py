# -*- coding: utf-8 -*-

# Cluster and compare incidence of activation and inhibition across species

import settings as S
import pandas as pd
import os
import numpy as np
import pdb
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns
plt.ion()
plt.close('all')

# Minimum number of interactions required to print data
minval = 3


ki = S.read_cache('inhibiting')
act = S.read_cache('activating')
tax = S.read_cache('TaxonomicData_temp')

# Drop entries without organism
ki = ki[pd.notnull(ki['Organism'])]
act = act[pd.notnull(act['Organism'])]

# Convert LigandID to string
ki['LigandID'] = ki['LigandID'].astype(str)
act['LigandID'] = act['LigandID'].astype(str)

# Drop null values
ki = ki[pd.notnull(ki['LigandID'])]
act = act[pd.notnull(act['LigandID'])]

# We don't want duplicate measurements of the same EC:LigandID in the same organism
ki.index = [':'.join( [ki.at[row,'EC_number'],ki.at[row,'LigandID'],ki.at[row,'Organism']] ) for row in ki.index]

act.index = [':'.join([act.at[row,'EC_number'], act.at[row,'LigandID'], act.at[row,'Organism'] ]) for row in act.index]

ki = ki[~ki.index.duplicated()]
act = act[~act.index.duplicated()]

# Make tables
print('Cross tabulating...')
kitab = pd.crosstab(ki.EC_number, ki.LigandID)
acttab = pd.crosstab(act.EC_number, act.LigandID)

# Drop indices where row or column sums equal zero
kitab = kitab.loc[(kitab.sum(axis=1) > minval), (kitab.sum(axis=0) >minval)]
acttab = acttab.loc[(acttab.sum(axis=1) > minval), (acttab.sum(axis=0) >minval)]

print('Writing to file...')
kitab.to_csv('../cache/inh_crosstab.csv')
acttab.to_csv('../cache/act_crosstab.csv')