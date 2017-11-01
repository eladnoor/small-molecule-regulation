#*- coding: utf-8 -*-
"""
Created on Sun Oct  9 17:37:42 2016

@author: noore
"""

from plot_figures_for_paper import FigurePlotter, abs_eps_x_v, CONDITIONS
import matplotlib.pyplot as plt
import settings
import numpy as np
import seaborn as sns

sns.set(style='white')

if __name__ == "__main__":
    fp = FigurePlotter()

    fp.plot_figS5()

    ki = fp.ki
    
    ki = ki[ki['growth condition'].isin(CONDITIONS)]
    
    #%%
    fig, axs = plt.subplots(1, 1, figsize=(3, 3))
    ax = axs
    s_range = np.logspace(-3, 3, 1000) # 10 uM - 100 mM

    interaction = 'pep:2.7.1.11'

    ax.plot(s_range, list(map(abs_eps_x_v, s_range)), color='k', alpha=0.5,
            zorder=1)
    ax.set_title('inhibitors', fontsize=12)

    ax.set_xscale('log')
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(-1e-2, 1+1e-2)
    
    met, ec = interaction.split(':')
    ax.set_title('Inhibition of %s on %s' % (met.upper(), ec))
    sat_vs_elast = ki.loc[ki['met:EC'] == interaction]
    sat_vs_elast['abs(elasticity)'] = sat_vs_elast['elasticity'].abs()
    sat_vs_elast['I_over_KI'] = sat_vs_elast['concentration']/sat_vs_elast['KI_Value']
    sat_vs_elast.plot.scatter(x='I_over_KI', y='abs(elasticity)', ax=ax, alpha=1,
                              zorder=2)
    for cond, row in sat_vs_elast.iterrows():
        ax.text(1.2*row['I_over_KI'], row['abs(elasticity)'],
                row['growth condition'], ha='left', va='center', fontsize=8)
    
    ax.set_xlabel('inhibitor saturation $I/K_I$', fontsize=12)
    ax.set_ylabel('$|\epsilon_v^x|$')
    settings.savefig(fig, 'figP1')
    
#    km_pivoted.index = km_pivoted.index.str.upper()
#    ki_pivoted.index = ki_pivoted.index.str.upper()
#
#    # keep and reorder only the conditions that were pre-selected
#    km_pivoted = km_pivoted.loc[:, CONDITIONS]
#    ki_pivoted = ki_pivoted.loc[:, CONDITIONS]
#
#    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(18, 30))
#    sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
#                cbar=False, vmin=-1, vmax=1, cmap=settings.HEATMAP_COLORMAP, fmt='.2f')
#    ax0.set_xticklabels(list(km_pivoted.columns), fontsize=12, rotation=90)
#    ax0.set_yticklabels(reversed(km_pivoted.index), rotation=0, fontsize=6)
#    ax0.set_title('substrates', fontsize=20)
#    ax0.set_xlabel('growth condition', fontsize=16)
#    ax0.set_ylabel('')
#
#    clb1 = matplotlib.colorbar.make_axes(ax1)
#    sns.heatmap(ki_pivoted, ax=ax1, mask=ki_pivoted.isnull(),
#                cbar=True, vmin=-1, vmax=1, annot=True, cmap=settings.HEATMAP_COLORMAP,
#                cbar_ax=clb1[0], fmt='.2f')
#    ax1.set_xticklabels(list(ki_pivoted.columns), fontsize=12, rotation=90)
#    ax1.set_title('inhibitors', fontsize=20)
#    ax1.set_yticklabels(reversed(ki_pivoted.index),
#                        rotation=0, fontsize=10)
#    ax1.set_xlabel('growth condition', fontsize=16)
#    ax1.set_ylabel('')
#    clb1[0].set_ylabel('elasticity', fontsize=16)
#
#    settings.savefig(fig, 'figS5')
#    km_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
#                                   'km_elasticity_full.csv'))
#    ki_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
#                                   'ki_elasticity_full.csv'))
