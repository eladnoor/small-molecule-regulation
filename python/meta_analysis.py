# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:11:08 2016

@author: noore

This script will count the number of interactions per EC-metabolite pair
across all organisms.
"""
import settings as S
import pandas as pd
import os
import numpy as np

ki = S.read_cache('ki')
act = S.read_cache('activating')

ki_merge = ki.groupby(['EC_number', 'Compound'])
ki_count = ki_merge.count()['Organism'].reset_index()

act_merge = act.groupby(['EC_number', 'Compound'])
act_count = act_merge.count()['Organism'].reset_index()
                    
# Write the merged files with some extra fun info
ki_idx = zip(ki_count['EC_number'],ki_count['Compound'])
ki_count['UniqueOrganisms_Ki'] = [','.join( np.unique( ki.ix[ki_merge.groups[ item ],'Organism']) ) for item in ki_idx]

act_idx = zip(act_count['EC_number'],act_count['Compound'])
act_count['UniqueOrganisms_Act'] = [','.join( np.unique( act.ix[act_merge.groups[ item ],'Organism']) ) for item in act_idx]


count_df = pd.merge(ki_count, act_count, on=['EC_number', 'Compound'],
                    how='outer', suffixes=('_inh', '_act')).fillna(0)                    
count_df.to_csv(os.path.join(S.RESULT_DIR, 'ec_met_organism_counts.csv'))

print "Enzymes which are inhibited by a compound in at least 40 organisms:"
print count_df[count_df['Organism_inh'] > 40]

print "Enzymes which are activated by a compound in at least 40 organisms:"
print count_df[count_df['Organism_act'] > 40]