#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:43:18 2017

@author: noore
"""

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import os
from settings import RESULT_DIR

rcParams['font.family'] = 'sans-serif'
rcParams['mathtext.sf'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'
fig, axs = plt.subplots(2, 2, figsize=(7, 7))

# first, plot the MM rate law (as a function of s)

Vmax = 1 # umol/min
Km = 1 # mM
s_range = np.linspace(0, 50, 100) # 10 uM - 100 mM
v = lambda s: Vmax * s / (Km + s)
eps = lambda s: 1 - s / (Km + s)

arrowprops = dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4)

x_low = 0.1
x_high = 45.0

fig.text(0.5, 0.95, 'Michaelis-Menten kinetics', fontsize=17, ha='center')
fig.text(0.5, 0.47, 'non-competitive inhibition', fontsize=17, ha='center')
###############################################################################
ax = axs[0, 0]
ax.plot(s_range, list(map(v, s_range)), '-')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(x_low, v(x_low)), xycoords='data',
            xytext=(50, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(x_high, v(x_high)), xycoords='data',
            xytext=(0, -30), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.set_title('rate law')
ax.annotate(r'$v = V^+ \, \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.5), xycoords='axes fraction', fontsize=14,
            ha='center', va='center')

###############################################################################
ax = axs[0, 1]
ax.plot(s_range, list(map(eps, s_range)), '-')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('elasticity $\epsilon_s^v$')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(x_low, eps(x_low)),
            xytext=(50, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(x_high, eps(x_high)),
            xytext=(0, 30), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| = 1 - \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.5), xycoords='axes fraction', fontsize=14,
            ha='center', va='center')
ax.set_title('substrate elasticity')

###############################################################################
v = lambda s: Vmax * (1 - s / (Km + s))
eps = lambda s: -s / (Km + s)

ax = axs[1, 0]
ax.plot(s_range, list(map(v, s_range)), '-')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(x_low, v(x_low)),
            xytext=(50, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(x_high, v(x_high)),
            xytext=(0, 30), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.set_title('rate law')
ax.annotate(r'$v = V^+ ( 1 - \frac{I}{K_I + I} ) $', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.5), xycoords='axes fraction', fontsize=14,
            ha='center', va='center')

###############################################################################
ax = axs[1, 1]
ax.plot(s_range, list(map(eps, s_range)), '-')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('elasticity $\epsilon_I^v$')

ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(x_low, eps(x_low)),
            xytext=(50, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(x_high, eps(x_high)),
            xytext=(0, 30), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$\epsilon_I^v = -\frac{I}{K_I + I}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.5), xycoords='axes fraction', fontsize=14,
            ha='center', va='center')
ax.set_title('inhibitor elasticity')

###############################################################################
fig.tight_layout(pad=4, h_pad=5, w_pad=1)
fig.savefig(os.path.join(RESULT_DIR, 'mca_linear.svg'))
fig.savefig(os.path.join(RESULT_DIR, 'mca_linear.pdf'))