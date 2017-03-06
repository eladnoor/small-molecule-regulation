#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:52:53 2017

@author: noore
"""

from cobra import Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.solvers import cplex_solver as solver
import settings as S
import os

def find_futile_cycle(cobra_model, fname):
    """
        try to locate EGCs by blocking all transport reactions and
        maximizing the objective of creating ATP from ADP + Pi
    """
    model = cobra_model.copy()

    model.reactions.ATPM.lower_bound = 0
    model.reactions.ATPM.upper_bound = 0

    # disable exchange reactions and ATM maintenance
    for r in model.reactions:
        if r.id[0:3] == 'EX_':
            r.lower_bound = 0
            r.upper_bound = 0

    # protons are sometimes not balanced well, so we ignore them
    model.reactions.EX_h_e.lower_bound = -1000
    model.reactions.EX_h_e.upper_bound = 1000
    model.reactions.EX_h2o_e.lower_bound = -1000
    model.reactions.EX_h2o_e.upper_bound = 1000

    # add a reaction that provides a steady flux of ATP
    atp_gen = Reaction("ATP_generator")
    atp_gen.add_metabolites({model.metabolites.atp_c: 1,
                             model.metabolites.h2o_c: 1,
                             model.metabolites.adp_c: -1,
                             model.metabolites.pi_c: -1,
                             model.metabolites.h_c: -1})
    atp_gen.lower_bound = 1
    atp_gen.upper_bound = 1
    model.add_reaction(atp_gen)

    last_solution_rids = []

    with open(fname, 'w') as fp:
        while True:
            # minimize the sum of all reactions
            for rid in last_solution_rids:
                model.reactions.get_by_id(rid).upper_bound = 0
            lp = solver.create_problem(model)
            for i, reaction in enumerate(model.reactions):
                solver.change_variable_objective(lp, i, 1)

            solver.solve_problem(lp, objective_sense='minimize')
            solution = solver.format_solution(lp, cobra_model)
            if solution.status != 'optimal':
                break
            fp.write("futile cycle with sum of fluxes = %g: \n" % (solution.f-1))
            for rid, flux in solution.x_dict.iteritems():
                if rid == 'ATP_generator':
                    continue
                if abs(flux) > 1e-9:
                    fp.write('%20s: %6.2f\n' % (rid, flux))
                    # disable flux in this reaction for the following iterations
                    model.reactions.get_by_id(rid).upper_bound = 0

###############################################################################

model_name = '../data/iJO1366'; obj = 'Ec_biomass_iJO1366_core_53p95M'
model = create_cobra_model_from_sbml_file(model_name + '.xml')
convert_to_irreversible(model)

# use MILP to find futile cycles
fname = os.path.join(S.RESULT_DIR, 'futile_cycles.txt')
find_futile_cycle(model, fname)