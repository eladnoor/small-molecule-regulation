#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:43:18 2017

@author: noore
"""

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import settings

rcParams['font.family'] = 'sans-serif'
rcParams['mathtext.sf'] = 'serif'
rcParams['mathtext.fontset'] = 'cm'
fig, axs = plt.subplots(2, 2, figsize=(8, 7))

# first, plot the MM rate law (as a function of s)

Vmax = 1 # umol/min
Km = 1 # mM
s_range = np.logspace(-3, 3, 100) # 10 uM - 100 mM
v_s = lambda s: Vmax * s / (Km + s)
eps_s_v = lambda s: 1 - s / (Km + s)
v_x = lambda s: Vmax * (1 - s / (Km + s))
eps_x_v = lambda s: -s / (Km + s)
abs_eps_x_v = lambda s: s / (Km + s)

x_low = 1e-2
x_high = 1e2

arrowprops = dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4)

fig.text(0.5, 0.95, 'Michaelis-Menten kinetics', fontsize=17, ha='center')
fig.text(0.5, 0.47, 'Non-competitive inhibition', fontsize=17, ha='center')

###############################################################################
ax = axs[0, 0]
ax.plot(s_range, map(v_s, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(x_low, v_s(x_low)),
            xytext=(60, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(x_high, v_s(x_high)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.set_title('rate law')
ax.annotate(r'$v = V^+ \, \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.1), xycoords='axes fraction', fontsize=14)

###############################################################################
ax = axs[0, 1]
ax.plot(s_range, map(eps_s_v, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlabel('substrate conc. $s$ [mM]')
ax.set_ylabel('elasticity $|\epsilon_s^v|$')

ax.annotate(r'$|\epsilon_s^v| \approx 1$', xy=(x_low, eps_s_v(x_low)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| \approx 0$', xy=(x_high, eps_s_v(x_high)),
            xytext=(0, 40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_s^v| = 1 - \frac{s}{K_M + s}$', color=(0.2, 0.4, 1.0),
            xy=(0.05, 0.1), xycoords='axes fraction', fontsize=14)
ax.set_title('substrate elasticity')

###############################################################################
ax = axs[1, 0]
ax.plot(s_range, map(v_x, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('rate $v$ [$\mu$mol/min]')

ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(x_low, v_x(x_low)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(x_high, v_x(x_high)),
            xytext=(-60, 0), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4))
ax.set_title('rate law')
ax.annotate(r'$v = V^+ ( 1 - \frac{I}{K_I + I} ) $', color=(0.2, 0.4, 1.0),
            xy=(0.05, 0.1), xycoords='axes fraction', fontsize=14)

###############################################################################
ax = axs[1, 1]
ax.plot(s_range, map(abs_eps_x_v, s_range), '-')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlabel('inhibitor conc. $I$ [mM]')
ax.set_ylabel('elasticity $|\epsilon_I^v|$')

ax.annotate(r'$|\epsilon_I^v| \approx 0$', xy=(x_low, abs_eps_x_v(x_low)),
            xytext=(0, 40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_I^v| \approx 1$', xy=(x_high, abs_eps_x_v(x_high)),
            xytext=(0, -40), textcoords='offset points', va='center', ha='center',
            fontsize=12,
            arrowprops=arrowprops)
ax.annotate(r'$|\epsilon_I^v| = \frac{I}{K_I + I}$', color=(0.2, 0.4, 1.0),
            xy=(0.5, 0.1), xycoords='axes fraction', fontsize=14)
ax.set_title('inhibitor elasticity')

###############################################################################
fig.tight_layout(pad=4, h_pad=5, w_pad=1)
settings.savefig(fig, 'mca_log', dpi=300)

###############################################################################
# make a small version of the figure to be placed as a subplot in the heatmap
# figures
#%%
fig, axs = plt.subplots(1, 2, figsize=(4, 2), sharey=True)
#axs[0].set_facecolor('#BBBBBB')
#axs[1].set_facecolor('#BBBBBB')

s_range = np.logspace(-3, 3, 1000) # 10 uM - 100 mM
eps = map(eps_s_v, s_range)
axs[0].plot([1e-3, 1e3], [0, 0], '--', color=(0.8, 0.8, 0.8))
axs[0].scatter(s_range, eps, c=eps, cmap=settings.HEATMAP_COLORMAP,
               edgecolor='none', s=15, vmin=-1, vmax=1)
eps = map(eps_x_v, s_range)
axs[1].plot([1e-3, 1e3], [0, 0], '--', color=(0.8, 0.8, 0.8))
axs[1].scatter(s_range, eps, c=eps, cmap=settings.HEATMAP_COLORMAP,
               edgecolor='none', s=15, vmin=-1, vmax=1)
axs[0].set_title('substrates')
axs[1].set_title('inhibitors')
axs[0].set_xlabel('substrate conc. $s$ [mM]')
axs[1].set_xlabel('inhibitor conc. $I$ [mM]')
axs[0].set_ylabel('elasticity')

axs[0].set_xscale('log')
axs[1].set_xscale('log')
axs[0].set_xlim(1e-3, 1e3)
axs[1].set_xlim(1e-3, 1e3)
axs[0].set_ylim(-1, 1)
settings.savefig(fig, 'elasticity_comparison', dpi=300)
