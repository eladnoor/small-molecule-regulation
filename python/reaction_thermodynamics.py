# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:37:50 2016

@author: noore

calculates dG'0 and dG'm for reaction from a Cobra model

"""

import settings

import os
import numpy as np
import pandas as pd
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T
from cobra.io.sbml import create_cobra_model_from_sbml_file

class reaction_thermodynamics(object):

    def __init__(self):
        reactions = reaction_thermodynamics.get_reactions_from_model()

        self.reactions = []
        self._not_balanced = []

        rstrings = []
        for r in reactions:
            k = r.kegg_reaction
            if k:
                if k.is_balanced() and not k.is_empty():
                    rstrings.append(k.write_formula())
                    self.reactions.append(r)
            else:
                self._not_balanced.append(r)
        self.Kmodel = KeggModel.from_formulas(rstrings)

        self.cc = ComponentContribution.init()
        self.pH = 7.5
        self.I = 0.2
        self.RT = R * default_T

    @staticmethod
    def get_reactions_from_model():
        cobra_model = create_cobra_model_from_sbml_file(settings.ECOLI_SBML_FNAME)
        metab_info = pd.DataFrame.from_csv(settings.ECOLI_MODEL_METABOLITES, sep='\t')        
        metab_info.dropna(inplace=True)
        for m in cobra_model.metabolites:
            try:
                m.CID = metab_info.kegg_id[m.id[:-2]]
            except KeyError:
                m.CID = None

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
                              columns=[r"dG'0", r"dG'0 std", r"dG'm", r"logRI"],
                              dtype=float)

        for i, r in enumerate(self.reactions):
            res_df.at[r.id.lower(), r"dG'0"] = dG0_prime[i, 0]
            res_df.at[r.id.lower(), r"dG'0 std"] = dG0_cov[i, 0]
            res_df.at[r.id.lower(), r"dG'm"] = dGm_prime[i, 0]
            res_df.at[r.id.lower(), r"logRI"] = logRI[i, 0]
        return res_df.round(2)
        
if __name__ == '__main__':
    Th = reaction_thermodynamics()
    df = Th.get_thermodynamics()
    df.to_csv(settings.ECOLI_THERMO_CACHE_FNAME)