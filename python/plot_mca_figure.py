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
s_range = np.logspace(-3, 3, 100) # 10 uM - 100 mM
v = lambda s: Vmax * s / (Km + s)
eps = lambda s: 1 - s / (Km + s)

arrowprops = dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4)
bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)

###############################################################################
ax = axs[0, 0]
ax.plot(s_range, map(v, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(2e-2, v(2e-2)),
            xytext=(60, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(1e2, v(1e2)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.set_title('Michaelis-Menten')
ax.annotate(r'$v = V^+ \cdot \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.1), xycoords='axes fraction', fontsize=14)

###############################################################################
ax = axs[0, 1]
ax.plot(s_range, map(eps, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('elasticity $\epsilon_s^v$')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(2e-2, eps(2e-2)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(1e2, eps(1e2)),
            xytext=(0, 40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| = 1 - \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.05, 0.1), xycoords='axes fraction', fontsize=14)
ax.set_title('substrate elasticity')

###############################################################################
v = lambda s: Vmax * (1 - s / (Km + s))
eps = lambda s: -s / (Km + s)

ax = axs[1, 0]
ax.plot(s_range, map(v, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(2e-2, v(2e-2)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(1e2, v(1e2)),
            xytext=(-60, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.set_title('non-competitive inhibition')
ax.annotate(r'$v = V^+ \cdot ( 1 - \frac{I}{K_I + I} ) $', color=(0.2, 0.4, 1.0),
            xy=(0.05, 0.1), xycoords='axes fraction', fontsize=14)

###############################################################################
ax = axs[1, 1]
ax.plot(s_range, map(eps, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('elasticity $\epsilon_I^v$')

ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(2e-2, eps(2e-2)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(1e2, eps(1e2)),
            xytext=(0, 40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$\epsilon_I^v = -\frac{I}{K_I + I}$', color=(0.2, 0.4, 1.0),
            xy=(0.05, 0.1), xycoords='axes fraction', fontsize=14)
ax.set_title('inhibitor elasticity')

###############################################################################
fig.tight_layout()
fig.savefig(os.path.join(RESULT_DIR, 'mca_log.pdf'))