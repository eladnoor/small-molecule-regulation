# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 23:42:54 2016

@author: eladn
"""

import settings
import os
import pandas as pd

ki = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_ki_bigg.csv'))
activators = pd.DataFrame.from_csv(os.path.join(settings.CACHE_DIR, 'ecoli_activating_bigg.csv'))

model_reactions = settings.get_reaction_table_from_xls()
bigg2ec = model_reactions.loc[:, ['Reaction Abbreviation', 'EC Number']]
bigg2ec.rename(columns={'Reaction Abbreviation': 'bigg.reaction'}, inplace=True)
bigg2ec = bigg2ec.loc[~bigg2ec['EC Number'].isnull()]

# change all reaction IDs to lower-case (apparently the standards have changed
# since the model was published, and cases are different now).

bigg2ec['bigg.reaction'] = bigg2ec['bigg.reaction'].apply(unicode.lower)