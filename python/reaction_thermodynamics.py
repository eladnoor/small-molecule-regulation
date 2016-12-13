# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:37:50 2016

@author: noore

calculates dG'0 and dG'm for reaction from a Cobra model

"""

import settings
import bigg
import kegg

import numpy as np
import pandas as pd
import json
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T
from cobra.io.sbml import create_cobra_model_from_sbml_file

class reaction_thermodynamics(object):

    def __init__(self):
        self.prepare_mappings()
        
        reactions = self.get_reactions_from_model()
        self.reactions = []
        self._not_balanced = []

        self.rstrings = []
        for r in reactions:
            k = r.kegg_reaction
            if k:
                if k.is_balanced() and not k.is_empty():
                    self.rstrings.append(k.write_formula())
                    self.reactions.append(r)
            else:
                self._not_balanced.append(r)
        self.Kmodel = KeggModel.from_formulas(self.rstrings)

        self.cc = ComponentContribution.init()
        self.pH = 7.3
        self.I = 0.25
        self.RT = R * default_T

    def prepare_mappings(self):
        b = bigg.BiGG()
        k = kegg.KEGG()
        
        met_df = k.kegg_df.join(b.metabolite_df, on='chebiID')
        self.bigg2kegg = met_df.groupby('bigg.metabolite')['KEGG_ID'].min().to_dict()
        self.bigg2ec = b.reaction_df.groupby('bigg.reaction')['EC_number'].min().to_dict()

        with open(settings.ECOLI_JSON_FNAME) as fp:
            ecoli_model = json.load(fp, encoding='UTF-8')
        
        self.bigg2subsystem = {}
        for r in ecoli_model['reactions']:
            rid = r['id'].lower()
            if 'subsystem' in r:
                self.bigg2subsystem[rid] = r['subsystem']
            else:
                self.bigg2subsystem[rid] = None

    def get_reactions_from_model(self):
        """
            Read all the reaction descriptions in the SBML file, and 
            map them to KEGG reaction (using another file of all
            E. coli metabolites and their BiGG and KEGG IDs)
        """
        cobra_model = create_cobra_model_from_sbml_file(settings.ECOLI_SBML_FNAME)
        for m in cobra_model.metabolites:
            m.CID = self.bigg2kegg.get(m.id[:-2], None)

        for r in cobra_model.reactions:
            CIDS = dict(zip(r.metabolites.keys(), map(lambda x: x.CID, r.metabolites.keys())))
            if None in CIDS.values():
                r.kegg_reaction = None
            else:
                sparse = {CIDS[m]:v for m,v in r.metabolites.iteritems()
                                            if CIDS[m]!='C00080'}        
                r.kegg_reaction = KeggReaction(sparse)
        
        return cobra_model.reactions

    def get_thermodynamics(self):
        '''
            Calculates the dG0 of a list of a reaction.
            Uses the component-contribution package (Noor et al) to estimate
            the standard Gibbs Free Energy of reactions based on 
            component contribution  approach and measured values (NIST and Alberty)

            Calculates the reversibility index (RI) of a reaction.
            The RI represent the change in concentrations of metabolites
            (from equal reaction reactants) that will make the reaction reversible.
            That is, the higher RI is, the more irreversible the reaction.
            A convenient threshold for reversibility is RI>=1000, that is a change of
            1000% in metabolite concentrations is required in ordeer to flip the
            reaction direction. 
        '''
        
        self.Kmodel.add_thermo(self.cc)
        
        temp = self.Kmodel.get_transformed_dG0(pH=self.pH, I=self.I, T=298.15)

        dG0_prime, dG0_cov = np.matrix(temp[0]), np.matrix(temp[1])
        
        conc = np.matrix(np.ones((len(self.Kmodel.cids), 1))) * 1e-3 # concentrations in M
        adj_mets = np.matrix(np.ones((len(self.Kmodel.cids), 1)))
        if 'C00001' in self.Kmodel.cids:
            j = self.Kmodel.cids.index('C00001')
            conc[j, 0] = 1
            adj_mets[j, 0] = 0

        dGm_prime = dG0_prime + self.RT * (self.Kmodel.S.T * np.log(conc))
        adj_count = np.abs(self.Kmodel.S.T) * adj_mets 
        inv_adj_count = np.matrix(np.diag( map(lambda x: 1.0/x, adj_count.flat) ))

        logRI = 2.0 * (inv_adj_count * dGm_prime) / (self.RT)
        
        res_df = pd.DataFrame(index=map(lambda r: r.id.lower(), self.reactions),
                              dtype=float)
        res_df['EC_number']     = map(self.bigg2ec.get, res_df.index)
        res_df['subsystem']     = map(self.bigg2subsystem.get, res_df.index)
        res_df['dG0_prime']     = dG0_prime
        res_df['dG0_prime_std'] = dG0_cov
        res_df['dGm_prime']     = dGm_prime
        res_df['logRI']         = logRI
        res_df['formula']       = self.rstrings
        res_df = res_df.round(2)
        return res_df
        
if __name__ == '__main__':
    Th = reaction_thermodynamics()
    df = Th.get_thermodynamics()
    df.to_csv(settings.ECOLI_THERMO_CACHE_FNAME)