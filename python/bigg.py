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
        self.ecoli_model, self.ecoli_S = settings.get_ecoli_json()
        self.metabolite_df = BiGG._get_metabolite_df()
        self.reaction_df = BiGG._get_reaction_df()

    @staticmethod
    def _get_metabolite_df():
        bigg2chebi = []
        with open(settings.BIGG_METABOLITE_FNAME, 'r') as fp:
            csv_reader = csv.reader(fp, delimiter='\t')
            next(csv_reader)
            for row in csv_reader:
                bigg_id = row[0]
                universal_bigg_id = row[1]
                database_links = json.loads(row[4])
                if 'CHEBI' in database_links:
                    for d in database_links['CHEBI']:
                        bigg2chebi.append(
                            (bigg_id, universal_bigg_id, d['id']))

        df = pd.DataFrame(bigg2chebi, columns=['bigg.metabolite with suffix',
                                               'bigg.metabolite', 'chebiID'])

        # keep only the first instance of each chebiID (to avoid
        # duplicating BRENDA data)
        return df[['bigg.metabolite', 'chebiID']].groupby('chebiID').first()

    @staticmethod
    def _get_reaction_df():
        bigg2ec = []

        # a few manual fixes for missing EC numbers in the model
        bigg2ec.append(('fhl', '1.1.99.33'))
        bigg2ec.append(('gdmane', '1.1.1.271'))

        # adding more EC numbers for citrate synthase
        # in the model there is only 2.3.3.3, which is citrate (Re)-synthase
        # we add 2.3.3.16 - citrate synthase (unknown stereospecificity)
        # and    2.3.3.1 - citrate (Si)-synthase
        # almost all data in BRENDA is for the unspecific EC number, therefore
        # it will be first
        bigg2ec.append(('cs', '2.3.3.16'))
        bigg2ec.append(('cs', '2.3.3.1'))

        # adding the real EC number for Pyruvate Dehydrogenase (PDH)
        # since BiGG only has the number for Dihydrolipoyllysine-residue
        # acetyltransferase which is a transient binding of the acyl
        # group to a modified lysine on the enzyme complex
        bigg2ec.append(('pdh', '1.2.4.1'))

        with open(settings.BIGG_REACTION_FNAME, 'r') as fp:
            csv_reader = csv.reader(fp, delimiter='\t')
            next(csv_reader)
            for row in csv_reader:
                bigg_id = row[0].lower()
                database_links = json.loads(row[3])
                if 'EC Number' in database_links:
                    for d in database_links['EC Number']:
                        bigg2ec.append((bigg_id, d['id']))

        df = pd.DataFrame(bigg2ec, columns=['bigg.reaction', 'EC_number'])

        return df

    def get_native_EC_numbers(self):
        model, _ = settings.get_ecoli_json()
        # make a list of all the EC numbers which are mapped to a BiGG reaction
        # in the E. coli model, which in turn is mapped to at least one E. coli
        # gene
        rids_with_genes = set()
        for d in self.ecoli_model['reactions']:
            rid = d['id']
            rule = d['gene_reaction_rule']
            if re.search('b[0-9]+', rule) is not None:
                rids_with_genes.add(rid.lower())

        # use the self.bigg object to convert these BiGG IDs to EC numbers
        bigg_reactions = self.reaction_df.groupby('EC_number').first()
        native_ec = set(bigg_reactions[bigg_reactions['bigg.reaction'].isin(
                                       rids_with_genes)].index)

        # adding Citrate Synthase
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

    def get_reaction_subsystems(self):
        bigg2subsystem = {}
        for r in self.ecoli_model['reactions']:
            rid = r['id'].lower()
            if 'subsystem' in r:
                bigg2subsystem[rid] = r['subsystem']
            else:
                bigg2subsystem[rid] = None
        bigg2subsystem = pd.DataFrame(columns=['bigg.subsystem'],
                                      data=bigg2subsystem.values(),
                                      index=bigg2subsystem.keys())
        return bigg2subsystem

if __name__ == '__main__':
    bigg = BiGG()
