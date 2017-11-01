# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:15:58 2016

@author: noore
"""
import settings
import pandas as pd
import os
import wesci
logger = wesci.Logger(script_file=__file__, log_file_prefix="./prefix")

logger.add_input_files({'regulation':
    os.path.join(settings.CACHE_DIR, 'regulation.csv')})

ORGANISM = 'Escherichia coli'

#%%
ki_df = settings.read_cache('regulation')
ki_df = ki_df[ki_df['Mode'] == '-']
ki_df.KI_Value.replace(-999, None, inplace=True)
ki_df = ki_df[ki_df['Organism'] == ORGANISM]
ki_df.drop(['Organism', 'Compound', 'LigandID', 'LigandName', 'Mode', 'Mechanism'],
           axis=1, inplace=True)
ki_df = ki_df.drop_duplicates()

ccm_df = pd.read_csv(settings.ECOLI_CCM_FNAME, index_col=None)
ccm_df.set_index('EC_number', inplace=True)

# select only Ki values that involve CCM enzymes
ccm_inh = ki_df.join(ccm_df, on='EC_number', how='inner')
ccm_inh['type'] = 'KI'
ccm_inh.sort_values(['EC_number', 'bigg.metabolite', 'type'], inplace=True)

# filter out mutated enzymes
ccm_inh = ccm_inh[ccm_inh['Commentary'].str.find('mutant') == -1]
ccm_inh = ccm_inh[ccm_inh['Commentary'].str.find('mutation') == -1]

# group by enzyme and metabolite names and select the smallest reported KI
ccm_ki = ccm_inh.groupby(['EC_number', 'bigg.metabolite']).min().reset_index()
ccm_ki = ccm_ki[~pd.isnull(ccm_ki['KI_Value'])]
ccm_ki.drop('type', axis=1, inplace=True)

#%%
output_fname = os.path.join(settings.RESULT_DIR, 'ecoli_ccm_inhibition.csv')
ccm_inh.to_csv(output_fname)
logger.add_output_files({'ecoli_ccm_inhibition': output_fname})

output_fname = os.path.join(settings.RESULT_DIR, 'ecoli_ccm_ki.csv')
ccm_ki.to_csv(output_fname)
logger.add_output_files({'ecoli_ccm_ki': output_fname})

logger.log()