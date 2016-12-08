# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 11:27:17 2016

@author: noore
"""

import os
import json
import settings
import pandas as pd
import numpy as np
from pulp import LpProblem, LpVariable, lpSum, LpMinimize, LpMaximize, LpStatusOptimal, solvers, value

class EcoliModel(object):
    
    PULP_SOLVER = None
    #PULP_SOLVER = solvers.CPLEX
    
    BIOMASS_REACTION = 'BIOMASS_Ec_iJO1366_WT_53p95M'
    #BIOMASS_REACTION = 'BIOMASS_Ec_iJO1366_core_53p95M'    
    
    def __init__(self):
        with open(settings.ECOLI_JSON_FNAME) as fp:
            self.model = json.load(fp, encoding='UTF-8')
            
        self.reactions = [d['id'] for d in self.model['reactions']]
        self.metabolites = [d['id'] for d in self.model['metabolites']]
        self.S = pd.DataFrame(index=self.metabolites, columns=self.reactions, dtype=float)
    
        for reaction in self.model['reactions']:
            r_column = self.S[reaction['id']]
            for met, coeff in reaction['metabolites'].iteritems():
                r_column[met] = float(coeff)
        self.S.fillna(0, inplace=True)
        
    def _CreateLinearProblem(self):
        """
            Generates a linear problem object (using PuLP) with all the mass
            balance constraints (Sv = 0), and with individual upper and lower
            bounds on reactions (as given by the model)
        """
        lp = LpProblem("FBA", LpMaximize)
        
        # create the flux variables
        fluxes = LpVariable.dicts("v", self.reactions,
                                  lowBound=-1000, upBound=1000)
                             
        for d in self.model['reactions']:
            if d['id'] not in self.reactions:
                print d['id'], "is not found in the list of reactions"
            fluxes[d['id']].lowBound = d['lower_bound']
            fluxes[d['id']].upBound = d['upper_bound']
    
        # add the mass balance constraints for all metabolites
        mass_balance_constraints = {}
        for met in self.metabolites:
            row = self.S.loc[met, :]
            mul = [row[i] * fluxes[self.reactions[i]] for i in row.nonzero()[0]]
            mass_balance_constraints[met] = (lpSum(mul) == 0)
            lp += mass_balance_constraints[met], "mass_balance_%s" % met
        
        return lp, fluxes, mass_balance_constraints
        
    def Solve(self):
        lp, fluxes, mass_balance_constraints = self._CreateLinearProblem()

        # set the objective to be the biomass function
        lp.setObjective(fluxes[EcoliModel.BIOMASS_REACTION])

        lp.solve(EcoliModel.PULP_SOLVER)
        if lp.status != LpStatusOptimal:
            raise solvers.PulpSolverError("cannot solve MDF primal")
        v_sol = pd.DataFrame(data=[fluxes[r].value() for r in self.reactions],
                             index=self.reactions, columns=['flux'])
        shadow_prices = pd.DataFrame(data=[mass_balance_constraints[m].pi for m in self.metabolites],
                                     index=self.metabolites, columns=['shadow price'])
        
        return lp.objective.value(), v_sol, shadow_prices

    def SolveForReaction(self, reaction, max_growth_yield=None):
        if max_growth_yield is None:
            max_gr, _, _ = self.Solve()
            max_growth_yield = 0.99*max_gr
        
        lp, fluxes, mass_balance_constraints = self._CreateLinearProblem()
        
        # set the objective to be the selected reaction
        lp.setObjective(fluxes[reaction])
        
        # set a biomass constraint
        lp += fluxes[EcoliModel.BIOMASS_REACTION] >= max_growth_yield

        lp.solve(EcoliModel.PULP_SOLVER)
        if lp.status != LpStatusOptimal:
            raise solvers.PulpSolverError("cannot solve MDF primal")
        v_all = pd.DataFrame(data=[fluxes[r].value() for r in self.reactions],
                             index=self.reactions, columns=['flux'])
        shadow_prices = pd.DataFrame(data=[mass_balance_constraints[m].pi for m in self.metabolites],
                                     index=self.metabolites, columns=['shadow price'])
        
        return lp.objective.value(), v_all, shadow_prices

