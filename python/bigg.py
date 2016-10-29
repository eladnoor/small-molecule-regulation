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
        self.ecoli_model, self.ecoli_S = settings.get_ecoli_json()

    @staticmethod
    def _get_metabolite_df():
        bigg2chebi = []
        with open(settings.BIGG_METABOLITE_FNAME, 'r') as fp:
            csv_reader = csv.reader(fp, delimiter='\t')
            csv_reader.next()
            for row in csv_reader:
                bigg_id = row[0]
                universal_bigg_id = row[1]
                database_links = json.loads(row[4])
                if 'CHEBI' in database_links:
                    for d in database_links['CHEBI']:
                        bigg2chebi.append((bigg_id, universal_bigg_id, d['id']))
        
        df = pd.DataFrame(bigg2chebi, columns=['bigg.metabolite with suffix',
                                               'bigg.metabolite', 'chebiID'])
    
        # keep only the first instance of each chebiID (to avoid
        # duplicating BRENDA data)
        return df[['bigg.metabolite', 'chebiID']].groupby('chebiID').first()
    
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

    def get_native_EC_numbers(self):
        model, _ = settings.get_ecoli_json()
        # make a list of all the EC numbers which are mapped to a BiGG reaction 
        # in the E. coli model, which in turn is mapped to at least one E. coli gene
        rids_with_genes = set()       
        for d in self.ecoli_model['reactions']:
            rid = d['id']
            rule = d['gene_reaction_rule']
            if re.search('b[0-9]+', rule) is not None:
                rids_with_genes.add(rid.lower())
        
        # use the self.bigg object to convert these BiGG IDs to EC numbers
        bigg_reactions = self.reaction_df
        native_ec = set(bigg_reactions[bigg_reactions['bigg.reaction'].isin(rids_with_genes)].index)
        native_ec.add('2.3.3.16')
        return native_ec
    
    def get_mets_in_cytosol(self):
        # make a list of all the metabolites that are cytoplasmic
        mets_in_cytosol = set()
        
        # load the BiGG model for E. coli in order to define the full scope
        # of metabolites and reactions in this organism
        for d in self.ecoli_model['metabolites']:
            met = d['id']
            if d['compartment'] == u'c':
                mets_in_cytosol.add(met.lower()[:-2])
        return mets_in_cytosol

if __name__ == '__main__':
    bigg = BiGG()
    