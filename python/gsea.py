#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 18:36:10 2017

@author: noore
"""

import settings
import pandas as pd
import os
import numpy as np
from scipy.stats import fisher_exact as fishertest
import matplotlib.pyplot as plt
from matplotlib import gridspec


def GSEA(feature, draw_figure=False):

    N = len(feature)  # total number of qualifying genes
    M = sum(feature)  # number of qualifying in-group genes
    if M == 0:
        return (1, [])

    PValues = []
    for i in xrange(N):
        A = sum(feature[:i])
        B = sum(~feature[:i])
        C = sum(feature[i+1:])
        D = sum(~feature[i+1:])
        oddsratio, p_value = fishertest([[A, B], [C, D]], 'greater');
        PValues.append(p_value)

    PValues = np.array(PValues)
    enrichment = min(PValues)
    I = np.argmin(PValues)

    if draw_figure:
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1],  wspace=0.0, hspace=0.0)

        ax0 = plt.subplot(gs[0])
        ax0.plot(range(N), -np.log10(PValues), '-')
        ax0.set_ylabel(r'$-\log_{10}$(p-value)')
        ax0.plot(I, -np.log10(enrichment), 'or')
        ax0.text(I - 0.02*N, -np.log10(enrichment),
                 'p-value = %.2f' % enrichment,
                 ha='right')
        ax0.set_xlim(0, N)
        ax0.set_ylim(0, -np.log10(enrichment)*1.2)
        ax0.get_xaxis().set_visible(False)

        ax1 = plt.subplot(gs[1])
        ax1.bar(range(N), feature, 0.4)
        ax1.set_ylim(0, 1)
        ax1.set_xlim(0, N)
        ax1.set_xlabel('Set #')
        ax1.get_yaxis().set_visible(False)

    return enrichment, I, fig

###############################################################################

reg_thermo_df = pd.DataFrame.from_csv(os.path.join(settings.RESULT_DIR,
                                                   'reg_thermo.csv'))

reg_thermo_df.sort_values(reg_thermo_df.columns[2], inplace=True)
feature = (reg_thermo_df['bigg.metabolite'] == 0)

enrichment, I, fig = GSEA(feature, True)
settings.savefig(fig, 'figS7')