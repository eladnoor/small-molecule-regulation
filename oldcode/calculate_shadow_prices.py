# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:46:31 2016

@author: eladn
"""

from ecoli_model import EcoliModel
import os
import settings as S
import pandas as pd

m = EcoliModel()
growth_yield, v_all, pi = m.Solve()

YIELD_FACTOR = 0.1
max_growth_yield = growth_yield * YIELD_FACTOR

result_df = pd.DataFrame(index=m.metabolites, columns=m.reactions)
for i, r in enumerate(m.reactions):
    v_r, v_all, pi = m.SolveForReaction(r, max_growth_yield)
    result_df[r] = pi
    print '\r(%3d of %3d) %40s' % (i+1, len(m.reactions), r)

result_df = result_df.round(3)
S.write_cache('shadow_prices_%g' % YIELD_FACTOR, result_df)