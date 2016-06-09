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

ki = S.read_cache('ki')
act = S.read_cache('activating')

ki_count = ki.groupby(['EC_number', 'Compound']).count()['Organism'].reset_index()
act_count = act.groupby(['EC_number', 'Compound']).count()['Organism'].reset_index()

count_df = pd.merge(ki_count, act_count, on=['EC_number', 'Compound'],
                    how='outer', suffixes=('_inh', '_act')).fillna(0)
                    
count_df.to_csv(os.path.join(S.RESULT_DIR, 'ec_met_organism_counts.csv'))

print "Enzymes which are inhibited by a compound in at least 40 organisms:"
print count_df[count_df['Organism_inh'] > 40]

print "Enzymes which are activated by a compound in at least 40 organisms:"
print count_df[count_df['Organism_act'] > 40]