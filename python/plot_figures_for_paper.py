#*- coding: utf-8 -*-
"""
Created on Sun Oct  9 17:37:42 2016

@author: noore
"""
from bigg import BiGG
from kegg import KEGG
import settings
import cache
import colorsys
import sys
from distutils.util import strtobool

import pandas as pd
import os
import json
import seaborn as sns
import numpy as np
from scipy.stats import gmean, ranksums
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
from matplotlib import rcParams
import pdb  # this is a reminder for Elad not to remove this pdb import
from topology import calculate_distances
import itertools

sns.set('paper', style='white')

ORGANISM = 'Escherichia coli'

STAT_TABLE_INDEX = ['all entries',
                    'keeping only E. coli data',
                    'filtering out data about mutated enzymes',
                    'keeping only data mapping to BiGG model',
                    'unique metabolite-enzyme pairs',
                    'unique metabolites',
                    'unique enzymes']

N_ACT_LABEL = 'Number of activating interactions'
N_INH_LABEL = 'Number of inhibiting interactions'

CONDITIONS = ['Glucose', 'Fructose', 'Galactose', 'Gluconate', 'Mannitol',
              'Sorbitol', 'Mannose', 'Glycerol', 'Pyruvate', 'Lactate',
              'Acetate', 'Succinate', 'glcNAc']

GENERAL_COLOR   = '#939598'
CCM_COLOR       = '#556B2f'

METABOLITE_COLOR = sns.color_palette('Set2')[3]
ACTIVATOR_COLOR = sns.color_palette('Set2')[0] # green
SUBSTRATE_COLOR = sns.color_palette(settings.HEATMAP_COLORMAP)[-1]
INHIBITOR_COLOR = sns.color_palette(settings.HEATMAP_COLORMAP)[0]
BOTH_COLOR      = sns.color_palette('Set2')[5]

# Michaelis-Menten
Vmax = 1 # umol/min
Km = 1 # mM
s_range = np.logspace(-3, 3, 100) # 10 uM - 100 mM
v_s = lambda s: Vmax * s / (Km + s)
eps_s_v = lambda s: 1 - s / (Km + s)
v_x = lambda s: Vmax * (1 - s / (Km + s))
eps_x_v = lambda s: -s / (Km + s)
abs_eps_x_v = lambda s: s / (Km + s)


class FigurePlotter(object):

    def __init__(self, rebuild_cache=False):
        self.stat_df = pd.DataFrame(index=STAT_TABLE_INDEX,
                                    columns=['km', 'KM_Value',
                                             'regulation', 'KI_Value'])
        self.kegg = KEGG()

        self.bigg = BiGG()
        self.native_mets = self.bigg.get_mets_in_cytosol()
        self.native_ECs = self.bigg.get_native_EC_numbers()

        self.get_data()

        _fname = os.path.join(settings.RESULT_DIR, 'ecoli_interactions.csv')
        self.regulation.to_csv(_fname)

    def get_kinetic_param(self, name, value_col, organism=ORGANISM):
        k = settings.read_cache(name)
        self.stat_df[name].iat[0] = k.shape[0]  # all entries
        self.stat_df[value_col].iat[0] = (k[value_col] > 0).sum()

        k = k[k['Organism'] == organism]
        self.stat_df[name].iat[1] = k.shape[0]  # filtered by organsim
        self.stat_df[value_col].iat[1] = (k[value_col] > 0).sum()

        k = k[(pd.isnull(k['Commentary'])) |
              ((k['Commentary'].str.find('mutant') == -1) &
               (k['Commentary'].str.find('mutation') == -1) &
               (k['Commentary'].str.find('variant') == -1) &
               (k['Commentary'].str.find('genetically engineered') == -1))]
        self.stat_df[name].iat[2] = k.shape[0]  # filtering mutants
        self.stat_df[value_col].iat[2] = (k[value_col] > 0).sum()

        # remove values with unmatched ligand
        k = k[pd.notnull(k['bigg.metabolite'])]
        k['bigg.metabolite'] = k['bigg.metabolite'].str.lower()

        return k

    def filter_non_native_interactions(self, k):
        k = k[k['bigg.metabolite'].isin(self.native_mets)]
        k = k[k['EC_number'].isin(self.native_ECs)]
        return k

    @staticmethod
    def calc_sat(k, value_col, conc_df, agg_type='gmean'):
        # filter missing Km or Ki values and -999 cases.
        k = k[k[value_col] > 0]

        # choose the minimum/median/gmean value among all repeats
        k = k.groupby(['EC_number', 'bigg.metabolite'])[value_col]
        if agg_type == 'minimum':
            k = k.min()
        elif agg_type == 'gmean':
            k = k.apply(gmean)
        elif agg_type == 'median':
            k = k.median()

        k = k.reset_index()
        # join data with measured concentrations
        k = k.join(conc_df, on='bigg.metabolite', how='inner')

        # melt table so each line will be a combination of EC,
        # substrate/inhibitor and growth condition
        k = pd.melt(k, id_vars=('EC_number', 'bigg.metabolite', value_col),
                    var_name='growth condition', value_name='concentration')

        k['saturation'] = k['concentration'] / (k['concentration'] +
                                                k[value_col])
        k['met:EC'] = k['bigg.metabolite'].str.cat(k['EC_number'], sep=':')
        return k

    @staticmethod
    def calc_agg_sat(k, agg_type='median', value_col='elasticity'):
        """
            calculates the [S]/K_S for all matching EC-metabolite pairs,
            in log2-fold-change.

            Input:
                K_df    - a DataFrame with three columns: EC_number,
                          bigg.metabolite, Value
                conc_df - a DataFrame with
        """
        k_grp = k.groupby(('bigg.metabolite', 'growth condition'))
        if agg_type == 'median':
            fc_med = k_grp.median()
        elif agg_type == 'gmean':
            fc_med = k_grp.agg(lambda x: gmean(list(x)))

        fc_med = fc_med[[value_col]].reset_index()
        fc_med = fc_med.pivot('bigg.metabolite', 'growth condition',
                              value_col)
        return fc_med.sort_index(axis=0)

    @staticmethod
    def get_subsystem_data():
        """
            Returns:
                - 1-to-many mapping BiGG Reaction IDs to cellular subsystems
                - many-to-many mapping of BiGG metabolites IDs to subsystems
        """
        with open(settings.ECOLI_JSON_FNAME) as fp:
            ecoli_model = json.load(fp, encoding='UTF-8')

        subsystem_data = []
        stoich_data = []
        for r in ecoli_model['reactions']:
            rid = r['id'].lower()
            if 'subsystem' in r:
                subsystem_data.append((rid, r['subsystem']))
            if 'metabolites' in r:
                for met, coeff in r['metabolites'].iteritems():
                    stoich_data.append((rid, met, coeff))

        reaction_subsystem_df = pd.DataFrame(
            subsystem_data,
            columns=('bigg.reaction', 'bigg.subsystem.reaction'))
        reaction_subsystem_df.set_index('bigg.reaction', inplace=True)

        stoich_df = pd.DataFrame(stoich_data,
                                 columns=('bigg.reaction',
                                          'bigg.metabolite', 'coeff'))
        # now associate every metabolite to subsystems by joining the two
        # tables
        metabolite_subsystem_df = stoich_df.join(
            reaction_subsystem_df, on='bigg.reaction')
        metabolite_subsystem_df.rename(
            columns={'bigg.subsystem.reaction': 'bigg.subsystem.metabolite'},
            inplace=True)
        metabolite_subsystem_df.drop('bigg.reaction', axis=1, inplace=True)
        metabolite_subsystem_df.drop('coeff', axis=1, inplace=True)
        metabolite_subsystem_df.drop_duplicates(inplace=True)

        # keep only cytoplasmic metabolites, and remove the suffix _c
        metabolite_subsystem_df = metabolite_subsystem_df[
            metabolite_subsystem_df['bigg.metabolite'].str[-2:] == '_c']

        # then remove the _c suffix and convert to lowercase
        metabolite_subsystem_df.loc[:, 'bigg.metabolite'] = \
            metabolite_subsystem_df['bigg.metabolite'].map(
            lambda s: s[0:-2].lower())

        return reaction_subsystem_df, metabolite_subsystem_df

    def get_data(self):
        _df = pd.DataFrame.from_csv(settings.ECOLI_METAB_FNAME)

        mean_data_cols = sum(_df.columns.str.findall('(.*\(mean\).*)'), [])
        std_data_cols = sum(_df.columns.str.findall('(.*\(std\).*)'), [])

        self.met_conc_mean = _df.loc[:, mean_data_cols]  # take only the data columns
        self.met_conc_std = _df.loc[:, std_data_cols]

        # remove the _c suffix from the compound names and convert to lowercase
        self.met_conc_mean.index = self.met_conc_mean.index.map(
            lambda s: s[0:-2].lower())
        self.met_conc_std.index = self.met_conc_mean.index.map(
            lambda s: s[0:-2].lower())

        colmap = dict(map(lambda x: (x, x[:-7]), self.met_conc_mean.columns))
        self.met_conc_mean.rename(columns=colmap, inplace=True)

        # for legacy reasons, also calculate the km and ki tables, without
        # filtering out the non-native EC reactions (in order to
        # make the full heatmap)
        km_raw_unfiltered = self.get_kinetic_param('km', 'KM_Value')
        self.km_unfiltered_ALL = km_raw_unfiltered

        self.km_unfiltered = FigurePlotter.calc_sat(
            km_raw_unfiltered, 'KM_Value', self.met_conc_mean)

        regulation_unfiltered = self.get_kinetic_param(
            'regulation', 'KI_Value')

        self.ki_unfiltered_ALL = regulation_unfiltered

        ki_raw_unfiltered = regulation_unfiltered[
            ~pd.isnull(regulation_unfiltered['KI_Value'])]

        self.ki_unfiltered = FigurePlotter.calc_sat(
            ki_raw_unfiltered, 'KI_Value', self.met_conc_mean)

        km_raw = self.filter_non_native_interactions(km_raw_unfiltered)

        self.regulation = self.filter_non_native_interactions(
            regulation_unfiltered)

        self.calc_unique_stats(km_raw, 'km', 'KM_Value')
        self.calc_unique_stats(self.regulation, 'regulation', 'KI_Value')

        # choose only one bigg.reaction for each EC number (arbitrarily)
        ec2bigg = self.bigg.reaction_df.groupby('EC_number').first()

        self.km = FigurePlotter.calc_sat(km_raw, 'KM_Value',
                                         self.met_conc_mean)
        self.km['elasticity'] = 1.0 - self.km['saturation']
        self.km = self.km.join(ec2bigg, on='EC_number', how='left')

        self.ki = FigurePlotter.calc_sat(
            self.regulation[~pd.isnull(self.regulation['KI_Value'])],
            'KI_Value', self.met_conc_mean)
        self.ki['elasticity'] = -self.ki['saturation']

        self.ki = self.ki.join(ec2bigg, on='EC_number', how='left')

        self.regulation = self.regulation.join(ec2bigg,
                                               on='EC_number', how='left')

        # write out SMRN prior to mapping to subsystems
        self.regulation.to_csv(os.path.join(settings.CACHE_DIR,
                                            'iJO1366_SMRN.csv'), index=False)

        self.reaction_subsystem_df, self.metabolite_subsystem_df = \
            FigurePlotter.get_subsystem_data()
        self.regulation = self.regulation.join(self.reaction_subsystem_df,
                                               on='bigg.reaction', how='left')
        self.regulation = pd.merge(self.regulation,
                                   self.metabolite_subsystem_df,
                                   on='bigg.metabolite', how='left')

        self.ki.to_csv(os.path.join(settings.RESULT_DIR,
                                    'ki_saturation_full.csv'))
        self.km.to_csv(os.path.join(settings.RESULT_DIR,
                                    'km_saturation_full.csv'))
        self.stat_df.drop('km', axis=1, inplace=True)
        self.stat_df.to_csv(os.path.join(settings.RESULT_DIR,
                                         'statistics.csv'))

    def calc_unique_stats(self, k, name, value_col):
        self.stat_df[name].iat[3] = k.shape[0]
        self.stat_df[name].iat[4] = \
            k.groupby(('bigg.metabolite', 'EC_number')).first().shape[0]
        self.stat_df[name].iat[5] = \
            k.groupby('bigg.metabolite').first().shape[0]
        self.stat_df[name].iat[6] = k.groupby('EC_number').first().shape[0]

        k_val = k[k[value_col] > 0]
        self.stat_df[value_col].iat[3] = k_val.shape[0]
        self.stat_df[value_col].iat[4] = \
            k_val.groupby(('bigg.metabolite', 'EC_number')).first().shape[0]
        self.stat_df[value_col].iat[5] = \
            k_val.groupby('bigg.metabolite').first().shape[0]
        self.stat_df[value_col].iat[6] = \
            k_val.groupby('EC_number').first().shape[0]

    def plot_fig4(self):
        """
            Panels a-b are for testing the  hypothesis that irreversible
            reactions are more likely to be regulated allosterically.
            Panels c-f show the difference between the distributions of
            substrate-enzyme interactions and regulator-enzyme interactions
            in terms of Km/Ki, saturation and elasticity.
        """
        fig, axs = plt.subplots(3, 2, figsize=(6, 9))

        # get the irreversibility constant (absolute) for all reactions
        # in the BiGG iJO1336 model
        thermo_df = pd.DataFrame.from_csv(settings.ECOLI_THERMO_CACHE_FNAME)

        # remove data about reactions with std=0 (i.e. known values)
        # and reactions with std > 20 (high uncertainty)
        thermo_df = thermo_df[(thermo_df["dG0_prime_std"] > 0) &
                              (thermo_df["dG0_prime_std"] < 20)]

        # select the median value of log(gamma) for each EC number
        # (in general, there should be only one value for each
        # EC number anyway)
        irr_index_l = r"$| log_{10}(\Gamma) |$"
        thermo_df[irr_index_l] = thermo_df['log10(RI)'].abs()
        thermo_df = thermo_df[~pd.isnull(thermo_df.EC_number)]

        # print the regulation table joined with the irreversibility values
        _temp_df = self.regulation.join(thermo_df[irr_index_l],
                                        on='EC_number')
        _temp_df.to_csv(os.path.join(settings.RESULT_DIR,
                                     'regulation_with_thermo.csv'))

        # group the thermo table by EC number and subsystem, while
        # taking the median value of the irreversibility index
        reg_thermo_df = thermo_df.groupby(['EC_number', 'subsystem'])
        reg_thermo_df = reg_thermo_df[irr_index_l].median().reset_index()

        # count how many unique interaction each EC number has
        # counting by metabolites (ignoring the modes)
        reg_count_df = self.regulation.groupby('EC_number')['bigg.metabolite'].nunique()
        reg_thermo_df = reg_thermo_df.join(reg_count_df, on='EC_number', how='left')
        reg_thermo_df.fillna(0, inplace=True)

        reg_thermo_df['num_regulators'] = ''
        reg_thermo_df.loc[reg_thermo_df['bigg.metabolite'] == 0, 'num_regulators'] = '0 regulators'
        reg_thermo_df.loc[reg_thermo_df['bigg.metabolite'].isin((1, 2)), 'num_regulators'] = '1-2 regulators'
        reg_thermo_df.loc[reg_thermo_df['bigg.metabolite'] > 2, 'num_regulators'] = '3+ regulators'

        reg_thermo_df['Regulation'] = ''
        reg_thermo_df.loc[reg_thermo_df['bigg.metabolite'] == 0, 'Regulation'] = 'not regulated'
        reg_thermo_df.loc[reg_thermo_df['bigg.metabolite'] > 0, 'Regulation'] = 'regulated'

        reg_thermo_df.to_csv(os.path.join(settings.RESULT_DIR, 'reg_thermo.csv'))

        ccm_thermo_df = reg_thermo_df[
            reg_thermo_df.subsystem.isin(settings.CCM_SUBSYSTEMS)]

        ccm_thermo_df.to_csv(os.path.join(settings.RESULT_DIR,
                             'CCM_thermodynamics.csv'))

        sns.set_palette('Set2', 8, 1)
        ax = axs[0, 0]
        FigurePlotter.comparative_cdf(x='num_regulators', y=irr_index_l,
                                      data=reg_thermo_df, ax=ax,
                                      title='all E. coli reactions')
        ax.set_xlim(0, 10)
        ax.plot([3, 3], [0, 1], 'k:', alpha=0.3, linewidth=1)

        ax = axs[0, 1]
        FigurePlotter.comparative_cdf(x='num_regulators', y=irr_index_l,
                                      data=ccm_thermo_df, ax=ax,
                                      title='only CCM reactions')
        ax.set_xlim(0, 10)
        ax.set_ylabel('')
        ax.plot([3, 3], [0, 1], 'k:', alpha=0.3, linewidth=1)

        # correlate irreversibility also with the number of references and
        # unique regulating metabolites

        num_refs = self.regulation.groupby(
             'bigg.reaction')['Literature'].nunique()
        ixrefs = num_refs.index.intersection(thermo_df.index)
        thermo_df['Num_Refs'] = 0
        thermo_df.loc[ixrefs, 'Num_Refs'] = num_refs.loc[ixrefs]

        num_regs = self.regulation.groupby(
            'bigg.reaction')['bigg.metabolite'].nunique()
        ixmets = num_regs.index.intersection(thermo_df.index)
        thermo_df['Num_Regs'] = 0
        thermo_df.loc[ixmets, 'Num_Regs'] = num_regs.loc[ixmets]
        thermo_df['is regulated'] = 'No'
        thermo_df.ix[thermo_df['Num_Regs'] > 0, 'is regulated'] = 'Yes'

        met_intersection = set(self.km['bigg.metabolite']).intersection(
            self.ki['bigg.metabolite'])
        km_inter = self.km[self.km['bigg.metabolite'].isin(met_intersection)]
        ki_inter = self.ki[self.ki['bigg.metabolite'].isin(met_intersection)]

        ax = axs[1, 0]

        concentrations = pd.melt(self.met_conc_mean)['value']
        concentrations = concentrations[~pd.isnull(concentrations)]

        sns.kdeplot(np.log10(concentrations), cumulative=False, ax=ax, bw=.25,
                    linewidth=2, color=METABOLITE_COLOR, legend=False)
        ax.set_xlim(-2.1, 2.1)
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
        ax.set_xlabel(r'$[S]$ (in mM)')
        ax.set_ylabel(r'Probability density')
        ax.set_title('Measured metabolite conc.')

        ax = axs[1, 1]
        km_values = km_inter.groupby(('met:EC')).first()['KM_Value']
        ki_values = ki_inter.groupby(('met:EC')).first()['KI_Value']
        sns.kdeplot(np.log10(km_values), cumulative=False,
                    ax=ax, bw=.25, color=SUBSTRATE_COLOR,
                    label='substrates (N = %d)' % km_values.shape[0],
                    linewidth=2)
        sns.kdeplot(np.log10(ki_values), cumulative=False,
                    ax=ax, bw=.25, color=INHIBITOR_COLOR,
                    label='inhibitors (N = %d)' % ki_values.shape[0],
                    linewidth=2)
        ax.set_xlim(-2.1, 2.7)
        ax.set_ylim(0, 0.7)
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
        ax.set_xlabel(r'$K_S$ (in mM)')
        ax.set_title(r'Measured $K_{\rm S}$ values')

        ranksum_res = ranksums(km_values, ki_values)
        ax.text(0.5, 0.8, '$p_{ranksum}$ < %.1g' % ranksum_res.pvalue,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)
        ax.legend(loc='upper right')

        # compare Km and Ki for the intersection of EC numbers

        ax = axs[2, 0]
        ki_saturation = ki_inter['saturation']
        ki_saturation = ki_saturation[~pd.isnull(ki_saturation)]
        km_saturation = km_inter['saturation']
        km_saturation = km_saturation[~pd.isnull(km_saturation)]
        sns.kdeplot(km_saturation, cumulative=False, ax=ax, bw=.1,
                    label='substrates (N = %d)' % km_saturation.shape[0],
                    linewidth=2, color=SUBSTRATE_COLOR)
        sns.kdeplot(ki_saturation, cumulative=False, ax=ax, bw=.1,
                    label='inhibitors (N = %d)' % ki_saturation.shape[0],
                    linewidth=2, color=INHIBITOR_COLOR)

        ax.grid(visible=False)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_xlabel(r'$\frac{[S]}{[S] + K_S}$')
        ax.set_ylabel(r'Probability density')
        ax.set_title(r'Saturation levels')
        ax.legend(loc='upper center')

        ranksum_res = ranksums(km_saturation, ki_saturation)
        ax.text(0.5, 0.8, '$p_{ranksum}$ < 10$^{%d}$' %
                np.ceil(np.log10(ranksum_res.pvalue)),
                horizontalalignment='center',
                verticalalignment='top',
                transform=ax.transAxes)

        ax = axs[2, 1]
        ki_elasticity = ki_inter['elasticity'].abs()
        ki_elasticity = ki_elasticity[~pd.isnull(ki_elasticity)]
        km_elasticity = km_inter['elasticity'].abs()
        km_elasticity = km_elasticity[~pd.isnull(km_elasticity)]
        sns.kdeplot(km_elasticity, cumulative=False, ax=ax, bw=.1,
                    label='substrates (N = %d)' % km_saturation.shape[0],
                    linewidth=2, color=SUBSTRATE_COLOR)
        sns.kdeplot(ki_elasticity, cumulative=False, ax=ax, bw=.1,
                    label='inhibitors (N = %d)' % ki_saturation.shape[0],
                    linewidth=2, color=INHIBITOR_COLOR)

        ax.grid(visible=False)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_xlabel(r'$|\epsilon_s^v|$')
        ax.set_title(r'Elasticities')
        ax.legend(loc='upper center')

        ranksum_res = ranksums(km_elasticity, ki_elasticity)
        ax.text(0.5, 0.8, '$p_{ranksum}$ < 10$^{%d}$' %
                np.ceil(np.log10(ranksum_res.pvalue)),
                horizontalalignment='center',
                verticalalignment='top',
                transform=ax.transAxes)

        for i, ax in enumerate(axs.flat):
            ax.annotate(chr(ord('a') + i), xy=(0.02, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=14)
        fig.tight_layout()

        settings.savefig(fig, 'fig4')


    def plot_fig5(self):
        """
            draw heat maps of the [S]/Ki and [S]/Km values across
            the 8 conditions
        """
        def count_values(k, value_col):
            _tmp = k.groupby(('bigg.metabolite', 'growth condition'))
            _tmp = _tmp.count().reset_index().groupby('bigg.metabolite').max()
            return _tmp[value_col].apply(int)

        fig = plt.figure(figsize=(10, 8))
        gs1 = GridSpec(1, 2)
        gs1.update(left=0.2, right=0.8, top=0.95, bottom=0.7, wspace=0.2)
        ax1 = plt.subplot(gs1[0, 0])
        ax2 = plt.subplot(gs1[0, 1])

        gs2 = GridSpec(1, 1)
        gs2.update(left=0.15, right=0.9, top=0.6, bottom=0.15, wspace=0.1)
        ax3 = plt.subplot(gs2[0])
        axs = [ax1, ax2, ax3]

        s_range = np.logspace(-3, 3, 1000) # 10 uM - 100 mM
        eps = map(eps_s_v, s_range)
        axs[0].plot([1e-3, 1e3], [0, 0], '--', color=(0.8, 0.8, 0.8))
        axs[0].scatter(s_range, eps, c=eps, cmap=settings.HEATMAP_COLORMAP,
                       edgecolor='none', s=15, vmin=-1, vmax=1)
        eps = map(eps_x_v, s_range)
        axs[1].plot([1e-3, 1e3], [0, 0], '--', color=(0.8, 0.8, 0.8))
        axs[1].scatter(s_range, eps, c=eps, cmap=settings.HEATMAP_COLORMAP,
                       edgecolor='none', s=15, vmin=-1, vmax=1)
        axs[0].set_title('substrates', fontsize=12)
        axs[1].set_title('inhibitors', fontsize=12)
        axs[0].set_xlabel('substrate conc. $s$ [mM]', fontsize=12)
        axs[1].set_xlabel('inhibitor conc. $I$ [mM]', fontsize=12)
        axs[0].set_ylabel('elasticity', fontsize=12)

        axs[0].set_xscale('log')
        axs[1].set_xscale('log')
        axs[0].set_xlim(1e-3, 1e3)
        axs[1].set_xlim(1e-3, 1e3)
        axs[0].set_ylim(-1, 1)
        axs[1].set_ylim(-1, 1)

        km_sat_agg = FigurePlotter.calc_agg_sat(self.km)
        ki_sat_agg = FigurePlotter.calc_agg_sat(self.ki)

        # keep and reorder only the conditions that were pre-selected
        km_sat_agg = km_sat_agg.loc[:, CONDITIONS]
        ki_sat_agg = ki_sat_agg.loc[:, CONDITIONS]

        # count how many K_M/K_I values we have for each metabolite
        # (i.e. how many different EC numbers)
        km_counts = count_values(self.km, 'KM_Value')
        ki_counts = count_values(self.ki, 'KI_Value')
        counts = pd.DataFrame([km_counts, ki_counts]).transpose()
        # make a dictionary mapping from the metabolite name to the same
        # name, followed by the counts (km, ki)
        index_mapping = {}
        for i, row in counts.iterrows():
            index_mapping[i] = '%s (%g,%g)' % (str(i).upper(), row['KM_Value'],
                                               row['KI_Value'])
        sat_joined = km_sat_agg.join(ki_sat_agg, how='inner',
                                     lsuffix='_sub', rsuffix='_inh')
        ind = sat_joined.mean(axis=1).sort_values(axis=0,
                                                  ascending=False).index
        sat_joined = sat_joined.reindex_axis(ind, axis=0)
        sat_joined.rename(index=index_mapping, inplace=True)

        sns.heatmap(sat_joined,
                    ax=axs[2], mask=sat_joined.isnull(), annot=True, fmt='.2f',
                    cbar=False, vmin=-1, vmax=1, cmap=settings.HEATMAP_COLORMAP,
                    annot_kws={'fontdict': {'fontsize': 8}})

        # change xtick labels back to the original strings
        # (without the suffixes) and increase the font size
        axs[2].set_xticklabels(list(km_sat_agg.columns) + list(ki_sat_agg.columns),
                               rotation=90, fontsize=12)

        # rotate the metabolite names back to horizontal, and increase
        # the font size
        axs[2].set_yticklabels(reversed(sat_joined.index), rotation=0, fontsize=10)

        axs[2].set_xlabel('growth condition', fontsize=10)
        axs[2].set_ylabel('')
        axs[2].set_title('as substrates' + ' '*50 + 'as inhibitors', fontsize=12)

        axs[2].axvline(sat_joined.shape[1]/2, 0, 1, color='r')
        settings.savefig(fig, 'fig5')

    def plot_figS5(self):
        def pivot_and_sort(k, sort_by='mean'):
            k_piv = k.pivot('met:EC', 'growth condition', 'elasticity')
            if sort_by == 'mean':
                ind = k_piv.mean(axis=1).sort_values(axis=0, ascending=True).index
            elif sort_by == 'index':
                ind = sorted(k_piv.index)
            k_piv = k_piv.reindex_axis(ind, axis=0)
            return k_piv

        km = self.km
        ki = self.ki
        km_pivoted = pivot_and_sort(km, sort_by='index')
        ki_pivoted = pivot_and_sort(ki, sort_by='index')
        km_pivoted.index = km_pivoted.index.str.upper()
        ki_pivoted.index = ki_pivoted.index.str.upper()

        # keep and reorder only the conditions that were pre-selected
        km_pivoted = km_pivoted.loc[:, CONDITIONS]
        ki_pivoted = ki_pivoted.loc[:, CONDITIONS]

        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(18, 30))
        sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
                    cbar=False, vmin=-1, vmax=1, cmap=settings.HEATMAP_COLORMAP, fmt='.2f')
        ax0.set_xticklabels(list(km_pivoted.columns), fontsize=12, rotation=90)
        ax0.set_yticklabels(reversed(km_pivoted.index), rotation=0, fontsize=6)
        ax0.set_title('substrates', fontsize=20)
        ax0.set_xlabel('growth condition', fontsize=16)
        ax0.set_ylabel('')

        clb1 = matplotlib.colorbar.make_axes(ax1)
        sns.heatmap(ki_pivoted, ax=ax1, mask=ki_pivoted.isnull(),
                    cbar=True, vmin=-1, vmax=1, annot=True, cmap=settings.HEATMAP_COLORMAP,
                    cbar_ax=clb1[0], fmt='.2f')
        ax1.set_xticklabels(list(ki_pivoted.columns), fontsize=12, rotation=90)
        ax1.set_title('inhibitors', fontsize=20)
        ax1.set_yticklabels(reversed(ki_pivoted.index),
                            rotation=0, fontsize=10)
        ax1.set_xlabel('growth condition', fontsize=16)
        ax1.set_ylabel('')
        clb1[0].set_ylabel('elasticity', fontsize=16)

        settings.savefig(fig, 'figS5')
        km_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
                                       'heatmap_km_saturation.csv'))
        ki_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
                                       'heatmap_ki_saturation.csv'))

    def plot_fig2ab(self):

        def venn3_sets(set_a, set_b, set_c, set_labels, ax):
            # order of values for Venn diagram:
            # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
            Abc = len(set_a.difference(set_b.union(set_c)))
            aBc = len(set_b.difference(set_a.union(set_c)))
            abC = len(set_c.difference(set_a.union(set_b)))
            ABc = len(set_a.intersection(set_b).difference(set_c))
            AbC = len(set_a.intersection(set_c).difference(set_b))
            aBC = len(set_b.intersection(set_c).difference(set_a))
            ABC = len(set_a.intersection(set_b).intersection(set_c))
            venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
                  set_labels=set_labels, ax=ax)

        print "found %d native interactions in %s" % \
            (self.regulation.shape[0], ORGANISM)

        ind_inh = self.regulation['Mode'] == '-'
        ind_act = self.regulation['Mode'] == '+'

        inh_met = set(self.regulation.loc[ind_inh, 'bigg.metabolite'])
        act_met = set(self.regulation.loc[ind_act, 'bigg.metabolite'])
        inh_ec = set(self.regulation.loc[ind_inh, 'EC_number'])
        act_ec = set(self.regulation.loc[ind_act, 'EC_number'])

        fig, axs = plt.subplots(1, 2, figsize=(7, 5))
        venn3_sets(inh_met, act_met, self.native_mets, ax=axs[0],
                   set_labels=('inhibitors', 'activators',
                               'E. coli metabolites (%d total)' %
                               len(self.native_mets)))
        venn3_sets(inh_ec, act_ec, self.native_ECs, ax=axs[1],
                   set_labels=('inhibited', 'activated',
                               'E. coli reactions (%d total)' %
                               len(self.native_ECs)))
        axs[0].annotate('a', xy=(0.02, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        axs[1].annotate('b', xy=(0.02, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        settings.savefig(fig, 'fig2ab')

        res = {'inhibitors': list(inh_met), 'activators': list(act_met),
               'all_metabolites': list(self.native_mets),
               'inhibited': list(inh_ec), 'activated': list(act_ec),
               'all_reactions': list(self.native_ECs)}
        _fname = os.path.join(settings.RESULT_DIR, 'venn_groups.json')
        with open(_fname, 'w') as fp:
            json.dump(res, fp, indent=4)

    def plot_fig2cd(self,highconfidence):
        """
            Draw 2D histograms of the number of activating and
            inhibiting reactions
            grouped by metabolite and grouped by reaction

            Highconfidence indicates that we should only use
            edges which have two or more literature references
        """

        highc_string = '_highconfidence' if highconfidence else  ''

        def plot_jointhist(data, xlabel, ylabel, xmax, ymax, highconfidence):
            """
                plot the histogram as a scatter plot with marginal histograms,
                and ensure empty bins are easily distinguishable from ones
                that have at least 1 hit.

                generally use xcrit = 12,ycrit = 5 unless using only high
                confidence interactions, then xcrit = 6, ycrit = 2
            """
            x = data[xlabel]
            y = data[ylabel]

            if highconfidence:
            	xcrit = 6
            	ycrit = 2
            else:
            	xcrit = 12
            	ycrit = 5

            # First, plot the scatter plot
            g = sns.JointGrid(x=x, y=y, size=4,
                              xlim=(-1, xmax+1), ylim=(-1, ymax+1))
            g = g.plot_joint(plt.scatter, alpha=0.2)
            plt.gcf()

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

            # annotate only unique points
            # annotate only unique points
            ann_df = data.drop_duplicates((xlabel, ylabel), keep=False)
            ann_df = ann_df[(ann_df[xlabel] > xcrit) | (ann_df[ylabel] > ycrit)]
            for i, row in ann_df.iterrows():
                    plt.annotate(i, xy=(row[xlabel], row[ylabel]),
                                 xytext=(row[xlabel]+1, row[ylabel]+1),
                                 ha='center', va='top', size=10,
                                 textcoords='data')

            # Next, plot the marginal histograms
            g = g.plot_marginals(sns.distplot, kde=False)

            return plt.gcf()


        # if using high confidence edges, then find those edges
        cols = ('bigg.reaction', 'bigg.metabolite')
        if highconfidence:
            reg = self.regulation
            reg = reg[reg['Source'] == 'BRENDA']
            reg['RefList'] = [item.split(',') if pd.notnull(item) else 0 for item in reg['Literature']]
            reg['ShortHand'] = reg['bigg.metabolite'].str.cat('-->' + reg['bigg.reaction'])
            reglit = reg.groupby(cols)

            highc = pd.DataFrame( columns = ['NumRef','Refs'] )
            for ii in reglit.groups.keys():
                ixs = reglit.groups[ii]
                tempref = reg.ix[ixs,'RefList']
                refs = np.unique(list(itertools.chain.from_iterable(tempref)))
                highc.ix[ii[1] + '-->' + ii[0],] = [len(refs),','.join(refs)]

            fig, axs = plt.subplots(1, 1, figsize=(6, 5))
            axs.hist( highc['NumRef'],bins = 0.5 + np.arange(0,highc['NumRef'].max()) )
            axs.set_xlabel('Number of Literature References')
            axs.set_ylabel('Number of SMRN edges')
            settings.savefig(fig, 'histogram_highconfidence_SMRN')
            highc.to_csv( os.path.join(settings.RESULT_DIR, 'histogram_highconfidence_SMRN.csv') )

        # join interaction table with bigg.reaction IDs
        # and keep only one copy of each reaction-metabolite pair
        if highconfidence:
            bigg_effectors = reg[reg['ShortHand'].isin(highc[highc['NumRef'] > 1].index)]
            bigg_effectors = bigg_effectors.groupby(cols).first().reset_index()
        else:
            bigg_effectors = self.regulation.groupby(cols).first().reset_index()

        bigg_effectors.drop('KI_Value', axis=1, inplace=True)

        # add columns for counting the positive and negative interactions
        bigg_effectors[N_ACT_LABEL] = 0
        bigg_effectors[N_INH_LABEL] = 0
        bigg_effectors.loc[bigg_effectors['Mode'] == '+', N_ACT_LABEL] = 1
        bigg_effectors.loc[bigg_effectors['Mode'] == '-', N_INH_LABEL] = 1

        grouped_by_met = bigg_effectors.groupby('bigg.metabolite').sum()
        grouped_by_rxn = bigg_effectors.groupby('bigg.reaction').sum()

        xmax = max(grouped_by_met[N_INH_LABEL].max(),
                   grouped_by_rxn[N_INH_LABEL].max())
        ymax = max(grouped_by_met[N_ACT_LABEL].max(),
                   grouped_by_rxn[N_ACT_LABEL].max())

        fig = plot_jointhist(grouped_by_met, N_INH_LABEL, N_ACT_LABEL,
                             xmax, ymax, highconfidence)
        fig.get_axes()[0].annotate('c', xy=(0.02, 0.98),
                     xycoords='axes fraction', ha='left', va='top',
                     size=20)
        settings.savefig(fig, 'fig2c' + highc_string)

        fig = plot_jointhist(grouped_by_rxn, N_INH_LABEL, N_ACT_LABEL,
                             xmax, ymax, highconfidence)
        fig.get_axes()[0].annotate('d', xy=(0.02, 0.98),
                     xycoords='axes fraction', ha='left', va='top',
                     size=20)
        settings.savefig(fig, 'fig2d' + highc_string)

    def plot_figS3(self, n_colors=20):

        def create_random_colormap(labels):
            np.random.seed(60)
            n = len(labels)
            h_values = np.random.permutation(n) / (1.0*n)
            s_values = 0.3 + 0.4 * (np.arange(n) % 2)
            v_values = 0.3 + 0.2 * (np.arange(n) % 3)

            s_values[0] = 0
            v_values[0] = 0.5
            colors = [colorsys.hsv_to_rgb(h, s, l) for h, s, l in
                      zip(h_values, s_values, v_values)]

            return dict(zip(labels, colors))

        cols = ['bigg.metabolite', 'bigg.reaction',
                'bigg.subsystem.reaction', 'Mode']
        reg_unique = self.regulation[cols].drop_duplicates()

        ylabel = 'bigg.metabolite'
        met_system_counter = reg_unique.groupby(
            (ylabel, 'bigg.subsystem.reaction', 'Mode')).count().reset_index()

        act_table = met_system_counter[met_system_counter['Mode'] == '+']
        act_table = act_table.pivot(index=ylabel,
                                    columns='bigg.subsystem.reaction',
                                    values='bigg.reaction')

        inh_table = met_system_counter[met_system_counter['Mode'] == '-']
        inh_table = inh_table.pivot(index=ylabel,
                                    columns='bigg.subsystem.reaction',
                                    values='bigg.reaction')

        act_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_activating.csv'))
        inh_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_inhibiting.csv'))

        # TODO: calculate the total number of reactions in each subsystem
        # and use it to normalize the number of interactions

        # for each subsystem in E. coli, draw a stacked bar plot
        # of the different metabolites regulating them
        # (giving unique colors only for the first 'n')
        # sort the subsystems by number of interactions
        inh_table = inh_table[inh_table.sum(0).sort_values().index]

        # lump some of the metabolites into more general groups
        lump_dict = {'divalent cations': ['CA2', 'CD2', 'CU2', 'MN2', 'MG2',
                                          'HG2', 'NI2', 'ZN2', 'FE2',
                                          'COBALT2'],
                     'monovalent cations': ['K', 'NA1', 'CL'],
                     'ATP/ADP/AMP': ['AMP', 'ADP', 'ATP'],
                     'GTP/GDP/GMP': ['GMP', 'GDP', 'GTP'],
                     'CTP/CDP/CMP': ['CMP', 'CDP', 'CTP'],
                     'UTP/UDP/UMP': ['UMP', 'UDP', 'UTP'],
                     'ITP/IDP/IMP': ['IMP', 'ITP'],
                     'Pi/PPi': ['PI', 'PPI'],
                     'NAD(P)(H)': ['NAD', 'NADH', 'NADP', 'NADPH'],
                     'amino acids': ['CYS__L', 'ALA__L', 'SER__L', 'GLU__L',
                                     'PRO__L', 'TYR__L', 'TRP__L', 'THR__L',
                                     'VAL__L', 'MET__L', 'ILE__L', 'HOM__L',
                                     'HIS__L', 'GLN__L', 'ARG__L', 'ASP__L'],
                     'cyanide': ['CYAN']
                     }
        inh_table_lumped = inh_table.transpose()
        inh_table_lumped.columns = map(str.upper, inh_table_lumped.columns)
        for k, v in lump_dict.iteritems():
            inh_table_lumped[k] = inh_table_lumped[v].sum(1)
            inh_table_lumped.drop(v, axis=1, inplace=True)

        # keep only the n metabolites with the most interactions
        ind = inh_table_lumped.sum(0).sort_values(ascending=False).index
        other = inh_table_lumped[ind[n_colors:]].sum(1).to_frame()
        counts_cutoff = pd.concat([other, inh_table_lumped[ind[0:n_colors]]],
                                  axis=1)
        counts_cutoff.rename(columns={0: 'other'}, inplace=True)

        colors = create_random_colormap(counts_cutoff.columns)
        colors['CYAN'] = (0.0, 0.8, 0.8)
        colors['UREA'] = (0.8, 0.8, 0.0)
        colors['Pi/PPi'] = (0.3, 0.8, 0.3)

        fig, ax = plt.subplots(1, 1, figsize=(6, 8))
        counts_cutoff.plot.barh(stacked=True, ax=ax,
                                color=map(colors.get, counts_cutoff.columns))
        ax.set_xlabel('# inhibitions')
        ax.set_ylabel('')
        fig.tight_layout()
        settings.savefig(fig, 'figS3', dpi=300)

    @staticmethod
    def cluster_matrix(X):
        import scipy.cluster.hierarchy as sch
        import scipy.spatial.distance as dist

        # X is a pandas dataframe
        X = X.fillna(0)
        Xd = dist.pdist(X, metric='euclidean')
        Xd2 = dist.squareform(Xd)
        links = sch.linkage(Xd2, method='ward', metric='euclidean')
        dendro = sch.dendrogram(links, no_plot=True)
        leaves = dendro['leaves']

        return leaves

    @staticmethod
    def comparative_cdf(x, y, data, ax, linewidth=2, title=None):
        xvals = sorted(data[x].unique())
        plt = sns.cubehelix_palette(len(xvals), start=2.3, rot=-.5,
                                    gamma=1.5, dark=0.4, light=0.7)
        for xval, c in zip(xvals, plt):
            d = data.loc[data[x] == xval, y]
            sns.kdeplot(d, cumulative=True, ax=ax, bw=.15,
                        linewidth=linewidth,
                        label=xval + ' (N = %d)' % d.shape[0],
                        color=c)
        ax.set_ylim(0, 1)
        ax.set_ylabel(r'Cumulative distribution')
        ax.set_xlabel(y)
        ax.plot([1e3, 1e3], [0, 1], '-', alpha=0.3)

        if title is None:
            title = x
        if len(xvals) == 2:
            ranksum_res = ranksums(data.loc[data[x] == xvals[0], y],
                                   data.loc[data[x] == xvals[1], y])
            ax.set_title(title + '\n$p_{ranksum}$ < %.1g' % ranksum_res.pvalue)
            return ranksum_res.pvalue
        else:
            ax.set_title(title)
            return None

    def plot_figS6(self, filter_using_model=False):

        from statsmodels.sandbox.stats.multicomp import multipletests as padjust
        import scipy.stats as st

        km = self.km if filter_using_model else self.km_unfiltered_ALL
        ki = self.ki if filter_using_model else self.ki_unfiltered_ALL

        # Get rid of non-positive entries
        km = km[km['KM_Value'] > 0]
        ki = ki[ki['KI_Value'] > 0]

        # Drop duplicates for multiple conditions
        km = km.drop_duplicates(subset = ['EC_number', 'bigg.metabolite','KM_Value'])
        ki = ki.drop_duplicates(subset = ['EC_number', 'bigg.metabolite','KI_Value'])

        res = pd.DataFrame()

        res['KI_Values'] = ki.groupby('bigg.metabolite')['KI_Value'].median()
        res['KI_Number'] = ki.groupby('bigg.metabolite')['EC_number'].nunique()
        res['KM_Values'] = km.groupby('bigg.metabolite')['KM_Value'].median()
        res['KM_Number'] = km.groupby('bigg.metabolite')['EC_number'].nunique()

        # Drop rows where we don't have data for both
        res = res.dropna()

        # Keep only metabolites with at least 2 measurements of each
        res = res[res['KI_Number'] > 1]
        res = res[res['KM_Number'] > 1]

        print "Found %d metabolites with more than 2 KMs and more than 2 KIs" % res.shape[0]

        res['PValue'] = np.nan

        # for each metabolite, if there is sufficient data, test
        for ii in res.index:
            kid = ki[ki['bigg.metabolite'] == ii]['KI_Value']
            kmd = km[km['bigg.metabolite'] == ii]['KM_Value']

            s,p = st.mannwhitneyu(kid, kmd)
            res.at[ii, 'PValue'] = p
            res['QValue'] = padjust(res['PValue'], method='fdr_bh')[1]
        res = res.sort_values('PValue')

        # Calculate log ratio of values
        res['Ratio'] = np.log2( res['KI_Values'] / res['KM_Values'] )
        res['size'] = 2 * (res['KI_Number'] + res['KM_Number'])

        # make a log-log scatter plots and mark the significant metabolites
        # in red
        c_not_significant = (0.8, 0.7, 1.0)
        c_significant = (0.9, 0.2, 0.2)

        fig, ax = plt.subplots(figsize=(5, 5),
                               subplot_kw={'xscale':'log', 'yscale':'log'})
        ax.scatter(res['KI_Values'], res['KM_Values'],
                   s=res['size'], color=c_not_significant)
        for ii, row in res.iterrows():
            if row['QValue'] < 0.1:
                ax.scatter(row['KI_Values'], row['KM_Values'],
                           color=c_significant, s=row['size'])
                ax.annotate(ii, xy=(row['KI_Values'], row['KM_Values']),
                            xytext=(5, 10), textcoords='offset points',
                            ha='left', va='center', color=c_significant)
        plt.xlabel('median $K_I$')
        plt.ylabel('median $K_M$')

        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3",
                             color='grey' )

        settings.savefig(fig, 'figS6')

    def plot_figS1(self):
        smrn = pd.read_csv(os.path.join(settings.CACHE_DIR,
                                        'iJO1366_SMRN.csv'), index_col=None)
        smrn_dist, all_distances = calculate_distances(smrn)
        smrn_merged = pd.merge(smrn, smrn_dist, on=['bigg.metabolite',
                                                    'bigg.reaction'])
        dist_mode_df = smrn_merged.groupby(('bigg.metabolite',
                                            'bigg.reaction', 'Mode')).first()
        dist_mode_df = dist_mode_df[['distance']].reset_index()
        inh_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '-', 'distance']
        act_distances = dist_mode_df.loc[dist_mode_df['Mode'] == '+', 'distance']

        Nmax = 9
        bins = range(Nmax+1)
        args = {'alpha': 1, 'normed': True, 'align': 'left', 'bins': bins,
                'linewidth': 0, 'rwidth': 0.8}
        fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=False)
        axs[0, 0].hist(all_distances, color=GENERAL_COLOR, **args)
        axs[1, 0].hist(smrn_dist['distance'], color=BOTH_COLOR, **args)
        axs[0, 1].hist(inh_distances, color=INHIBITOR_COLOR , **args)
        axs[1, 1].hist(act_distances, color=ACTIVATOR_COLOR, **args)

        for i in range(2):
            axs[i, 0].set_ylabel('Fraction of metabolite-enzyme pairs')
            axs[1, i].set_xlabel('Distance in # reactions between '
                                 'metabolite and enzyme')
        for i, ax in enumerate(axs.flat):
            #ax.annotate(chr(ord('a') + i), xy=(0.02, 0.98),
            #            xycoords='axes fraction', ha='left', va='top',
            #            size=20)
            ax.set_xticks(np.arange(Nmax+1))
            ax.set_xlim(-1, Nmax+1)

        axs[0, 0].annotate('all iJO1366 pairs', xy=(0.9, 0.9),
                        xycoords='axes fraction', ha='right', va='top',
                        size=14)
        axs[1, 0].annotate('all SMRN pairs', xy=(0.9, 0.9),
                        xycoords='axes fraction', ha='right', va='top',
                        size=14)
        axs[0, 1].annotate('only inhibition', xy=(0.9, 0.9),
                        xycoords='axes fraction', ha='right', va='top',
                        size=14)
        axs[1, 1].annotate('only activation', xy=(0.9, 0.9),
                        xycoords='axes fraction', ha='right', va='top',
                        size=14)

        settings.savefig(fig, 'figS1')
        smrn_dist.to_csv(os.path.join(settings.CACHE_DIR,
                                      'iJO1366_SMRN_dist.csv'), index=False)

    def plot_figS2(self):
        smrn = pd.read_csv(os.path.join(settings.CACHE_DIR,
                                        'iJO1366_SMRN.csv'), index_col=None)
        #rxn_hist = smrn.join(self.reaction_subsystem_df, on='bigg.reaction')
        #rxn_hist = self.reaction_subsystem_df.copy()
        rxn_hist = smrn.groupby('bigg.reaction')['bigg.metabolite'].nunique()
        rxn_hist = rxn_hist.reset_index().join(self.reaction_subsystem_df,
                                               on='bigg.reaction')
        #rxn_hist['no. regulators'].fillna(0, inplace=True)
        rxn_hist_ccm = rxn_hist[rxn_hist['bigg.subsystem.reaction'].isin(settings.CCM_SUBSYSTEMS)]

        _mets = self.metabolite_subsystem_df
        _mets = _mets[_mets['bigg.subsystem.metabolite'].isin(settings.CCM_SUBSYSTEMS)]
        ccm_mets = set(_mets['bigg.metabolite'])
        # metabolites are not associated to unique subsystems, but we will chose
        # one arbitrarily
        #met_hist = pd.DataFrame(index=all_mets)
        met_hist = smrn.groupby('bigg.metabolite')['bigg.reaction'].nunique().reset_index()
        met_hist_ccm = met_hist[met_hist['bigg.metabolite'].isin(ccm_mets)]

        Nmax = 15
        bins = range(1, Nmax+1)
        args = {'alpha': 1, 'normed': False, 'align': 'left', 'bins': bins,
                'linewidth': 0, 'rwidth': 0.8}
        fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=False)
        axs[0, 0].hist(rxn_hist['bigg.metabolite'], color=GENERAL_COLOR, **args)
        axs[1, 0].hist(met_hist['bigg.reaction'], color=GENERAL_COLOR, **args)
        axs[0, 1].hist(rxn_hist_ccm['bigg.metabolite'], color=CCM_COLOR, **args)
        axs[1, 1].hist(met_hist_ccm['bigg.reaction'], color=CCM_COLOR, **args)

        axs[0, 0].set_ylabel('No. of reactions')
        axs[1, 0].set_ylabel('No. of metabolites')
        for i in range(2):
            axs[1, i].set_xlabel('No. of interactions')
        for i, ax in enumerate(axs.flat):
            #ax.annotate(chr(ord('a') + i), xy=(0.02, 0.98),
            #            xycoords='axes fraction', ha='left', va='top',
            #            size=20)
            ax.set_xticks(bins)
            ax.set_xlim(0, Nmax+1)

        axs[0, 0].annotate('all regulated reactions (N = %d)' % rxn_hist.shape[0],
                           xy=(0.9, 0.9), size=10,
                           xycoords='axes fraction', ha='right', va='top')
        axs[1, 0].annotate('all regulating metabolites (N = %d)' % met_hist.shape[0],
                           xy=(0.9, 0.9), size=10,
                           xycoords='axes fraction', ha='right', va='top')
        axs[0, 1].annotate('only CCM reactions (N = %d)' % rxn_hist_ccm.shape[0],
                           xy=(0.9, 0.9), size=10,
                           xycoords='axes fraction', ha='right', va='top')
        axs[1, 1].annotate('only CCM metabolites (N = %d)' % met_hist_ccm.shape[0],
                           xy=(0.9, 0.9), size=10,
                           xycoords='axes fraction', ha='right', va='top')
        settings.savefig(fig, 'figS2')

    def plot_figS4(self):
        rcParams['font.family'] = 'sans-serif'
        rcParams['mathtext.sf'] = 'serif'
        rcParams['mathtext.fontset'] = 'cm'
        fig, axs = plt.subplots(2, 2, figsize=(8, 7))

        # first, plot the MM rate law (as a function of s)

        x_low = 1e-2
        x_high = 1e2

        arrowprops = dict(facecolor='black', shrink=0.01, width=1.5, headwidth=4)

        fig.text(0.5, 0.95, 'Michaelis-Menten kinetics', fontsize=17, ha='center')
        fig.text(0.5, 0.47, 'Non-competitive inhibition', fontsize=17, ha='center')

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

        fig.tight_layout(pad=4, h_pad=5, w_pad=1)
        settings.savefig(fig, 'figS4', dpi=300)

    def print_ccm_table(self):
        ccm_df = pd.DataFrame.from_csv(settings.ECOLI_CCM_FNAME,
                                       index_col=None)
        ccm_df.set_index('EC_number', inplace=True)

        # select only entries that involve CCM enzymes
        ccm_ki = self.regulation.join(ccm_df, on='EC_number', how='inner')
        ccm_ki = ccm_ki[~pd.isnull(ccm_ki['KI_Value'])]
        ccm_interactions = self.regulation.join(ccm_df,
                                                on='EC_number', how='inner')
        ccm_interactions['KI_Value'] = None

        ccm_concat = pd.concat([ccm_ki, ccm_interactions])
        ccm_concat.sort_values(
            ['EC_number', 'bigg.metabolite', 'Mode'], inplace=True)

        ccm_concat.to_csv(os.path.join(settings.RESULT_DIR, 'ccm_data.csv'))
        return ccm_concat

def user_yes_no_query(question, default=False):
    if default:
        sys.stdout.write('%s? [(yes)/no] ' % question)
    else:
        sys.stdout.write('%s? [yes/(no)] ' % question)

    while True:
        try:
            ri = raw_input().lower()
            if not ri:
                return default
            return strtobool(ri)
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')

###############################################################################
if __name__ == "__main__":
    plt.close('all')

    rebuild_cache = user_yes_no_query('Rebuild cache files', default=False)
    if rebuild_cache:
        cache.rebuild_cache()

    fp = FigurePlotter()

#     fp.plot_fig2ab()
#    fp.plot_fig2cd(highconfidence = True)
    fp.plot_fig4()
#     fp.plot_fig5()
#
#     fp.plot_figS1()
#     fp.plot_figS2()
#     fp.plot_figS3()
#     fp.plot_figS4()
#     fp.plot_figS5()
#     fp.plot_figS6()
#
#     fp.print_ccm_table()

    plt.close('all')
