# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:13:37 2016

@author: noore
"""

import settings
import pandas as pd
import json
import csv

def get_metabolite_df():
    bigg2chebi = []
    with open(settings.BIGG_METABOLITE_FNAME, 'r') as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        csv_reader.next()
        for row in csv_reader:
            bigg_id = row[0]
            database_links = json.loads(row[4])
            if 'CHEBI' in database_links:
                for d in database_links['CHEBI']:
                    bigg2chebi.append((bigg_id, d['id'][6:]))
    
    df = pd.DataFrame(bigg2chebi, columns=['bigg.metabolite', 'ChEBI'])

    # keep only metabolites with a _c suffix (i.e. cytoplasmic)
    keep_idx = df['bigg.metabolite'].apply(lambda s: s[-2:] == '_c')
    return df[keep_idx].groupby('ChEBI').first()

def get_reaction_df():
    bigg2ec = []
    with open(settings.BIGG_REACTION_FNAME, 'r') as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        csv_reader.next()
        for row in csv_reader:
            bigg_id = row[0].lower()
            database_links = json.loads(row[3])
            if 'EC Number' in database_links:
                for d in database_links['EC Number']:
                    bigg2ec.append((bigg_id, d['id']))
    
    df = pd.DataFrame(bigg2ec, columns=['bigg.reaction', 'EC_number'])

    # keep only metabolites with a _c suffix (i.e. cytoplasmic)
    return df.groupby('EC_number').first()
