# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:13:37 2016

@author: noore
"""

import settings
import pandas as pd
import json
import csv
import re

class BiGG(object):
    
    def __init__(self):
        
        self.metabolite_df = BiGG._get_metabolite_df()
        self.reaction_df = BiGG._get_reaction_df()
        
    @staticmethod
    def _get_metabolite_df():
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
    
    @staticmethod
    def _get_reaction_df():
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
    
        return df.groupby('EC_number').first()

if __name__ == '__main__':
    bigg = BiGG()
    
    model, S = settings.get_ecoli_json()    
    rids_with_genes = []        
    for d in model['reactions']:
        rid = d['id']
        rule = d['gene_reaction_rule']
        if re.findall('b[0-9]+', rule) != []:
            rids_with_genes.append(rid.lower())
            
    native_ec = set(bigg.reaction_df[bigg.reaction_df['bigg.reaction'].isin(rids_with_genes)].index)
    
    
    mets_in_cytoplasm = set()       
    for d in model['metabolites']:
        met = d['id']
        if d['compartment'] == u'c':
            mets_in_cytoplasm.add(met.lower())
    
