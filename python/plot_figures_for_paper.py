# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 17:37:42 2016

@author: noore
"""
from bigg import BiGG
from kegg import KEGG
import settings
import map_ligands

import pandas as pd
import os, json
import seaborn as sns
import numpy as np
from scipy.stats import gmean, ranksums
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import matplotlib
sns.set('paper', style='white')

ORGANISM = 'Escherichia coli'
SAT_FORMULA_S = r'$[S]/\left([S] + K_S\right)$'
SAT_FORMULA_M = r'$[S]/\left([S] + K_M\right)$'
SAT_FORMULA_I = r'$[S]/\left([S] + K_I\right)$'

STAT_TABLE_INDEX = ['all entries',
                    'keeping only E. coli data',
                    'filtering out data about mutated enzymes',
                    'keeping only data mapping to BiGG model',
                    'unique metabolite-enzyme pairs',
                    'unique metabolites',
                    'unique enzymes']

class FigurePlotter(object):

    def __init__(self, rebuild_cache=False):
        self.stat_df = pd.DataFrame(index=STAT_TABLE_INDEX,
                                    columns=['km', 'KM_Value', 'regulation', 'KI_Value'])
        self.kegg = KEGG()

        self.bigg = BiGG()
        self.native_mets = self.bigg.get_mets_in_cytosol()
        self.native_ECs = self.bigg.get_native_EC_numbers()

        if rebuild_cache:
            map_ligands.rebuild_cache()

        self.get_data()

    def get_kinetic_param(self, name, value_col, organism=ORGANISM):
        k = settings.read_cache(name)
        self.stat_df[name].iat[0] = k.shape[0]  # all entries
        self.stat_df[value_col].iat[0] = (k[value_col] > 0).sum()

        k = k[k['Organism'] == organism]
        self.stat_df[name].iat[1] = k.shape[0]  # filtered by organsim
        self.stat_df[value_col].iat[1] = (k[value_col] > 0).sum()

        k = k[(pd.isnull(k['Commentary'])) |
              ((k['Commentary'].str.find('mutant') == -1) &
               (k['Commentary'].str.find('mutation') == -1))]
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

        k['saturation'] = k['concentration'] / (k['concentration'] + k[value_col])
        k['met:EC'] = k['bigg.metabolite'].str.cat(k['EC_number'], sep=':')
        return k

    @staticmethod
    def calc_agg_sat(k, agg_type='gmean'):
        """
            calculates the [S]/K_S for all matching EC-metabolite pairs,
            in log2-fold-change.

            Input:
                K_df    - a DataFrame with three columns: EC_number, bigg.metabolite, Value
                conc_df - a DataFrame with
        """
        if agg_type == 'median':
            fc_med = k.groupby(('bigg.metabolite', 'growth condition')).median()
        elif agg_type == 'gmean':
            fc_med = k.groupby(('bigg.metabolite', 'growth condition')).agg(lambda x: gmean(list(x)))

        fc_med = fc_med[['saturation']].reset_index()
        fc_med = fc_med.pivot('bigg.metabolite', 'growth condition', 'saturation')
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

        reaction_subsystem_df = pd.DataFrame(subsystem_data,
                                             columns=('bigg.reaction', 'bigg.subsystem.reaction'))
        reaction_subsystem_df.set_index('bigg.reaction', inplace=True)

        stoich_df = pd.DataFrame(stoich_data,
                                 columns=('bigg.reaction', 'bigg.metabolite', 'coeff'))

        # now associate every metabolite to subsystems by joining the two tables

        metabolite_subsystem_df = stoich_df.join(reaction_subsystem_df, on='bigg.reaction')
        metabolite_subsystem_df.rename(columns={'bigg.subsystem.reaction': 'bigg.subsystem.metabolite'}, inplace=True)
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
        self.met_conc_mean = _df.iloc[:, 1:9]  # take only the data columns
        self.met_conc_std = _df.iloc[:, 10:]

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
        self.km_unfiltered = FigurePlotter.calc_sat(
            km_raw_unfiltered, 'KM_Value', self.met_conc_mean)

        regulation_unfiltered = self.get_kinetic_param(
            'regulation', 'KI_Value')

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

        self.km = FigurePlotter.calc_sat(km_raw, 'KM_Value', self.met_conc_mean)
        self.km = self.km.join(ec2bigg, on='EC_number', how='left')

        self.ki = FigurePlotter.calc_sat(self.regulation[~pd.isnull(self.regulation['KI_Value'])],
                                         'KI_Value', self.met_conc_mean)

        self.ki = self.ki.join(ec2bigg, on='EC_number', how='left')

        self.regulation = self.regulation.join(ec2bigg, on='EC_number', how='left')
        
        # write out SMRN prior to mapping to subsystems
        self.regulation.to_csv(os.path.join(settings.RESULT_DIR, 'iJO1366_SMRN.csv'))

        reaction_subsystem_df, metabolite_subsystem_df = FigurePlotter.get_subsystem_data()
        self.regulation = self.regulation.join(reaction_subsystem_df,
                                               on='bigg.reaction', how='left')
        self.regulation = pd.merge(self.regulation, metabolite_subsystem_df,
                                   on='bigg.metabolite', how='left')

        self.ki.to_csv(os.path.join(settings.RESULT_DIR, 'ki_saturation_full.csv'))
        self.km.to_csv(os.path.join(settings.RESULT_DIR, 'km_saturation_full.csv'))
        self.stat_df.drop('km', axis=1, inplace=True)
        self.stat_df.to_csv(os.path.join(settings.RESULT_DIR, 'statistics.csv'))
        

    def calc_unique_stats(self, k, name, value_col):
        self.stat_df[name].iat[3] = k.shape[0]
        self.stat_df[name].iat[4] = k.groupby(('bigg.metabolite', 'EC_number')).first().shape[0]
        self.stat_df[name].iat[5] = k.groupby('bigg.metabolite').first().shape[0]
        self.stat_df[name].iat[6] = k.groupby('EC_number').first().shape[0]

        k_val = k[k[value_col] > 0]
        self.stat_df[value_col].iat[3] = k_val.shape[0]
        self.stat_df[value_col].iat[4] = k_val.groupby(('bigg.metabolite', 'EC_number')).first().shape[0]
        self.stat_df[value_col].iat[5] = k_val.groupby('bigg.metabolite').first().shape[0]
        self.stat_df[value_col].iat[6] = k_val.groupby('EC_number').first().shape[0]

    def draw_agg_heatmaps(self, agg_type='median'):
        """
            draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
        """
        km_sat_agg = FigurePlotter.calc_agg_sat(self.km, agg_type)
        ki_sat_agg = FigurePlotter.calc_agg_sat(self.ki, agg_type)
        sat_joined = km_sat_agg.join(ki_sat_agg, how='inner', lsuffix='_sub', rsuffix='_inh')
        sat_joined = sat_joined.reindex_axis(sat_joined.mean(axis=1).sort_values(axis=0, ascending=False).index, axis=0)

        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        clb = matplotlib.colorbar.make_axes(ax)

        sns.heatmap(sat_joined, ax=ax, mask=sat_joined.isnull(), annot=True,
                    cbar=True, vmin=0, vmax=1, cmap='viridis', cbar_ax=clb[0],
                    annot_kws={'fontdict': {'fontsize': 12}})

        # change xtick labels back to the original strings (without the suffixes) and increase the font size
        ax.set_xticklabels(list(km_sat_agg.columns) + list(ki_sat_agg.columns), rotation=90, fontsize=12)

        # rotate the metabolite names back to horizontal, and increase the font size
        ax.set_yticklabels(reversed(sat_joined.index), rotation=0, fontsize=12)

        ax.set_xlabel('growth condition', fontsize=16)
        ax.set_ylabel('')
        ax.set_title('substrates' + ' '*30 + 'inhibitors', fontsize=20)
        clb[0].set_ylabel('%s saturation over all reactions' % agg_type, fontsize=16)
        clb[0].set_yticklabels(np.linspace(0.0, 1.0, 6), fontsize=12)

        ax.axvline(sat_joined.shape[1]/2, 0, 1, color='r')

        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation_%s.svg' % agg_type))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation_%s.png' % agg_type), dpi=300)

    def draw_full_heapmats(self, filter_using_model=True):
        if filter_using_model:
            km_pivoted = self.km.pivot('met:EC', 'growth condition', 'saturation')
            ki_pivoted = self.ki.pivot('met:EC', 'growth condition', 'saturation')
        else:
            km_pivoted = self.km_unfiltered.pivot('met:EC', 'growth condition', 'saturation')
            ki_pivoted = self.ki_unfiltered.pivot('met:EC', 'growth condition', 'saturation')

        km_pivoted = km_pivoted.reindex_axis(km_pivoted.mean(axis=1).sort_values(axis=0, ascending=True).index, axis=0)
        ki_pivoted = ki_pivoted.reindex_axis(ki_pivoted.mean(axis=1).sort_values(axis=0, ascending=True).index, axis=0)

        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 30))
        sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
                    cbar=False, vmin=0, vmax=1, cmap='viridis')
        ax0.set_xticklabels(list(km_pivoted.columns), fontsize=12, rotation=90)
        ax0.set_yticklabels(reversed(km_pivoted.index), rotation=0, fontsize=6)
        ax0.set_title('substrates', fontsize=20)
        ax0.set_xlabel('growth condition', fontsize=16)
        ax0.set_ylabel('')

        clb1 = matplotlib.colorbar.make_axes(ax1)
        sns.heatmap(ki_pivoted, ax=ax1, mask=ki_pivoted.isnull(),
                    cbar=True, vmin=0, vmax=1, annot=True, cmap='viridis', cbar_ax=clb1[0])
        ax1.set_xticklabels(list(ki_pivoted.columns), fontsize=12, rotation=90)
        ax1.set_title('inhibitors', fontsize=20)
        ax1.set_yticklabels(reversed(ki_pivoted.index), rotation=0, fontsize=10)
        ax1.set_xlabel('growth condition', fontsize=16)
        ax1.set_ylabel('')
        clb1[0].set_ylabel('saturation', fontsize=16)
        clb1[0].set_yticklabels(np.linspace(0.0, 1.0, 6), fontsize=12)

        if filter_using_model:
            fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation.svg'))
            fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation.png'), dpi=200)
            km_pivoted.to_csv(os.path.join(settings.RESULT_DIR, 'heatmap_km_saturation.csv'))
            ki_pivoted.to_csv(os.path.join(settings.RESULT_DIR, 'heatmap_ki_saturation.csv'))
        else:
            fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation_unfiltered.png'), dpi=100)

    def draw_venn_diagrams(self):

        def venn3_sets(set_a, set_b, set_c, set_labels, ax):
            # order of values for Venn diagram: (Abc, aBc, ABc, abC, AbC, aBC, ABC)
            Abc = len(set_a.difference(set_b.union(set_c)))
            aBc = len(set_b.difference(set_a.union(set_c)))
            abC = len(set_c.difference(set_a.union(set_b)))
            ABc = len(set_a.intersection(set_b).difference(set_c))
            AbC = len(set_a.intersection(set_c).difference(set_b))
            aBC = len(set_b.intersection(set_c).difference(set_a))
            ABC = len(set_a.intersection(set_b).intersection(set_c))
            venn3(subsets = (Abc, aBc, ABc, abC, AbC, aBC, ABC),
                  set_labels=set_labels, ax=ax)

        print "found %d native interactions in %s" % (self.regulation.shape[0], ORGANISM)

        ind_inh = self.regulation['Mode']=='-'
        ind_act = self.regulation['Mode']=='+'

        inh_met = set(self.regulation.loc[ind_inh, 'bigg.metabolite'])
        act_met = set(self.regulation.loc[ind_act, 'bigg.metabolite'])
        inh_ec = set(self.regulation.loc[ind_inh, 'EC_number'])
        act_ec = set(self.regulation.loc[ind_act, 'EC_number'])

        fig, axs = plt.subplots(1, 2, figsize=(7, 5))
        venn3_sets(inh_met, act_met, self.native_mets,
                   set_labels=('inhibitors', 'activators', 'E. coli metabolites'), ax=axs[0])
        venn3_sets(inh_ec, act_ec, self.native_ECs,
                   set_labels=('inhibited', 'activated', 'E. coli reactions'), ax=axs[1])
        axs[0].annotate('a', xy=(0.02, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        axs[1].annotate('b', xy=(0.02, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        fig.savefig(os.path.join(settings.RESULT_DIR, 'venn.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'venn.png'), dpi=600)

        res = {'inhibitors': list(inh_met), 'activators': list(act_met),
               'all_metabolites': list(self.native_mets),
               'inhibited': list(inh_ec), 'activated': list(act_ec),
               'all_reactions': list(self.native_ECs)}
        with open(os.path.join(settings.RESULT_DIR, 'venn_groups.json'), 'w') as fp:
            json.dump(res, fp, indent=4)
        self.regulation.to_csv(open(os.path.join(settings.RESULT_DIR, 'ecoli_interactions.csv'), 'w'))

    def draw_cdf_plots(self, linewidth=2):
        """
            Compare the CDFs of the two fold-change types (for Ki and Km)
        """

        met_color ='#fedf08' # dandelion
        ki_color = '#fe86a4' # rosa
        km_color = '#3eaf76' # dark seafoam green

        fig, axs = plt.subplots(1, 3, figsize=(7.5, 3), sharey=True)

        met_intersection = set(self.km['bigg.metabolite']).intersection(self.ki['bigg.metabolite'])
        km_inter = self.km[self.km['bigg.metabolite'].isin(met_intersection)]
        ki_inter = self.ki[self.ki['bigg.metabolite'].isin(met_intersection)]

        ax = axs[0]

        concentrations = pd.melt(self.met_conc_mean)['value']
        concentrations = concentrations[~pd.isnull(concentrations)]

        #concentrations.hist(cumulative=True, normed=1, bins=1000, histtype='step', ax=ax, linewidth=1, color='orange')
        sns.kdeplot(np.log10(concentrations), cumulative=True, ax=ax, bw=.15,
                    linewidth=linewidth, color=met_color, legend=False)
        ax.set_xlim(-2.1, 2.1)
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
        ax.set_ylim(0, 1)
        ax.set_xlabel(r'$[S]$ (in mM)')
        ax.set_ylabel(r'Cumulative distribution')
        ax.set_title('Measured metabolite conc.')

        ax = axs[1]
        sns.kdeplot(-np.log10(km_inter['KM_Value']), cumulative=True, ax=ax, bw=.15,
                    label='substrates $(K_M)$', linewidth=linewidth, color=km_color)
        sns.kdeplot(-np.log10(ki_inter['KI_Value']), cumulative=True, ax=ax, bw=.15,
                    label='inhibitors $(K_I)$', linewidth=linewidth, color=ki_color)
        ax.set_xlim(-2.1, 2.7)
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
        ax.set_ylim(0, 1)
        ax.set_xlabel(r'$K_S$ (in mM)')
        ax.set_title(r'Measured $K_{\rm S}$ values')

        ranksum_res = ranksums(km_inter['KM_Value'], ki_inter['KI_Value'])
        ax.text(0.5, 0.1, '$p_{ranksum}$ < %.1g' % ranksum_res.pvalue,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)
        ax.legend(loc='upper left')

        # compare Km and Ki for the intersection of EC numbers

        ax = axs[2]
        ki_saturation = ki_inter['saturation']
        ki_saturation = ki_saturation[~pd.isnull(ki_saturation)]
        km_saturation = km_inter['saturation']
        km_saturation = km_saturation[~pd.isnull(km_saturation)]
        sns.kdeplot(km_saturation, cumulative=True, ax=ax, bw=.01,
                    label='substrates', linewidth=linewidth, color=km_color)
        sns.kdeplot(ki_saturation, cumulative=True, ax=ax, bw=.01,
                    label='inhibitors', linewidth=linewidth, color=ki_color)
        ax.grid(visible=False)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0, 1)
        ax.set_xlabel(SAT_FORMULA_S)
        ax.set_title(r'Saturation levels')
        ax.legend(loc='upper left')

        ranksum_res = ranksums(km_saturation, ki_saturation)
        ax.text(0.5, 0.1, '$p_{ranksum}$ < %.1g' % ranksum_res.pvalue,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)
        fig.tight_layout()

        fig.savefig(os.path.join(settings.RESULT_DIR, 'saturation_histogram.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'saturation_histogram.png'), dpi=600)

    def draw_2D_histograms(self, tax2use='kingdom', minsize=10):
        """
            Draw 2D histograms of the number of activating and inhibiting reactions
            grouped by metabolite and grouped by reaction
        """

        def plot_2d_hist(data, xlabel, ylabel, ax, xmax, ymax):
            """
                plot the histogram as a heatmap, and make sure
                empty bins are easily distinguishable from ones that have
                at least 1 hit.
            """
            x = data[xlabel]
            y = data[ylabel]
            bins = (np.arange(0, xmax+2), np.arange(0, ymax+2))
            H, _, _ = np.histogram2d(x, y, bins=bins)
            # in a heatmap, the matrix columns become rows and vice versa
            # so we use H.T instead of H
            sns.heatmap(np.log2(H.T), mask=(H.T==0), ax=ax, annot=False,
                        cbar=True, vmin=0, vmax=8, cmap='viridis')
            ax.invert_yaxis()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_xticks(0.5+np.arange(0, xmax+1, 5))
            ax.set_xticklabels(np.arange(0, xmax+1, 5), rotation=0)

            cax = plt.gcf().axes[-1]
            cax.set_yticks(np.arange(0, 9, 1))
            cax.set_yticklabels(2**np.arange(0, 9, 1))

            for i in data.index:
                x_i = data.at[i, xlabel]
                y_i = data.at[i, ylabel]
                if (x_i > 12 or y_i > 5) and (H[x_i, y_i] == 1):
                    #print i, x_i, y_i, H[x_i, y_i]
                    ax.annotate(i, xy=(x_i, ymax-y_i), xytext=(x_i+1, ymax-y_i-1),
                                ha='center', va='top', size=6,
                                textcoords='data')

        # join interaction table with bigg.reaction IDs
        # and keep only one copy of each reaction-metabolite pair
        cols = ('bigg.reaction', 'bigg.metabolite')
        bigg_effectors = self.regulation.groupby(cols).first().reset_index()

        # add columns for counting the positive and negative interactions
        n_act_label = 'Number of activating interactions'
        n_inh_label = 'Number of inhibiting interactions'

        bigg_effectors[n_act_label] = 0
        bigg_effectors[n_inh_label] = 0
        bigg_effectors.loc[bigg_effectors['Mode'] == '+', n_act_label] = 1
        bigg_effectors.loc[bigg_effectors['Mode'] == '-', n_inh_label] = 1

        grouped_by_met = bigg_effectors.groupby('bigg.metabolite').sum()
        grouped_by_rxn = bigg_effectors.groupby('bigg.reaction').sum()

        # plot the data as 2D histograms

        fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=False)
        axs[0].annotate('grouped by metabolite', xy=(0.5, 0.98),
                        xycoords='axes fraction', ha='center',
                        va='top', size=12)
        axs[1].annotate('grouped by reaction', xy=(0.5, 0.98),
                        xycoords='axes fraction', ha='center',
                        va='top', size=12)
        axs[0].annotate('c', xy=(0.98, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        axs[1].annotate('d', xy=(0.98, 0.98),
                        xycoords='axes fraction', ha='left', va='top',
                        size=20)
        xmax = max(grouped_by_met[n_inh_label].max(), grouped_by_rxn[n_inh_label].max())
        ymax = max(grouped_by_met[n_act_label].max(), grouped_by_rxn[n_act_label].max())

        plot_2d_hist(grouped_by_met, n_inh_label, n_act_label, axs[0], xmax, ymax)
        axs[0].set_xlabel('')
        plot_2d_hist(grouped_by_rxn, n_inh_label, n_act_label, axs[1], xmax, ymax)

        fig.savefig(os.path.join(settings.RESULT_DIR, '2D_histograms.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, '2D_histograms.png'), dpi=600)

    def print_ccm_table(self):
        ccm_df = pd.DataFrame.from_csv(settings.ECOLI_CCM_FNAME, index_col=None)
        ccm_df.set_index('EC_number', inplace=True)

        # select only entries that involve CCM enzymes
        ccm_ki = self.regulation.join(ccm_df, on='EC_number', how='inner')
        ccm_ki = ccm_ki[~pd.isnull(ccm_ki['KI_Value'])]
        ccm_interactions = self.regulation.join(ccm_df, on='EC_number', how='inner')
        ccm_interactions['KI_Value'] = None

        ccm_concat = pd.concat([ccm_ki, ccm_interactions])
        ccm_concat.sort_values(['EC_number', 'bigg.metabolite', 'Mode'], inplace=True)

        ccm_concat.to_csv(os.path.join(settings.RESULT_DIR, 'ccm_data.csv'))
        return ccm_concat

    def draw_pathway_histogram(self):
        pathways = set(self.regulation['bigg.subsystem.reaction'])
        pathways.update(self.regulation['bigg.subsystem.metabolite'])
        pathways = sorted(pathways)
        if np.nan in pathways:
            pathways.remove(np.nan)

        cols = ['bigg.metabolite', 'bigg.subsystem.metabolite', 'bigg.reaction', 'bigg.subsystem.reaction', 'Mode']
        reg_unique = self.regulation[cols].drop_duplicates()

        ylabel = 'bigg.subsystem.metabolite'
        met_system_counter = reg_unique.groupby(
            (ylabel, 'bigg.subsystem.reaction', 'Mode')).count().reset_index()

        act_table = met_system_counter[met_system_counter['Mode'] == '+']
        act_table = act_table.pivot(index=ylabel,
                                    columns='bigg.subsystem.reaction',
                                    values='bigg.reaction').fillna(0)

        inh_table = met_system_counter[met_system_counter['Mode'] == '-']
        inh_table = inh_table.pivot(index=ylabel,
                                    columns='bigg.subsystem.reaction',
                                    values='bigg.reaction').fillna(0)

        N = len(pathways)
        hist_mat_act = np.zeros((N, N))
        for i in xrange(N):
            if pathways[i] in act_table.index:
                for j in xrange(N):
                    if pathways[j] in act_table.columns:
                        hist_mat_act[i, j] = act_table.at[pathways[i], pathways[j]]

        N = len(pathways)
        hist_mat_inh = np.zeros((N, N))
        for i in xrange(N):
            if pathways[i] in inh_table.index:
                for j in xrange(N):
                    if pathways[j] in inh_table.columns:
                        hist_mat_inh[i, j] = inh_table.at[pathways[i], pathways[j]]

        vmax = max(hist_mat_act.max(), hist_mat_inh.max())

        fig, axs = plt.subplots(1, 2, figsize=(20, 9))

        ax = axs[0]
        sns.heatmap(hist_mat_act, ax=ax, annot=False, cbar=False, cmap='viridis',
                    xticklabels=pathways, yticklabels=pathways, vmin=0, vmax=vmax)
        ax.set_title('activation')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('activated enzyme subsystem')

        ax = axs[1]
        sns.heatmap(hist_mat_inh, ax=ax, annot=False, cbar=True, cmap='viridis',
                    xticklabels=pathways, yticklabels=False, vmin=0, vmax=vmax)
        ax.set_title('inhibition')
        ax.set_xlabel('inhibited enzyme subsystem')

        fig.tight_layout(pad=4)
        fig.savefig(os.path.join(settings.RESULT_DIR, 'pathway_histograms.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'pathway_histograms.png'), dpi=300)
        act_table.to_csv(os.path.join(settings.RESULT_DIR, 'pathway_histograms_activating.csv'))
        inh_table.to_csv(os.path.join(settings.RESULT_DIR, 'pathway_histograms_inhibiting.csv'))

    def draw_pathway_met_histogram(self):
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

        fig, axs = plt.subplots(1, 2, figsize=(20, 30), sharey=True)

        ax = axs[0]
        sns.heatmap(act_table, ax=ax, annot=False, cbar=True, cmap='viridis',
                    vmin=1, vmax=12)
        ax.set_title('activation')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('activated enzyme subsystem')

        ax = axs[1]
        sns.heatmap(inh_table, ax=ax, annot=False, cbar=True, cmap='viridis',
                    vmin=1, vmax=12)
        ax.set_title('inhibition')
        ax.set_xlabel('inhibited enzyme subsystem')

        fig.tight_layout(pad=4)
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms.svg'))
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms.png'), dpi=300)
        act_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_activating.csv'))
        inh_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_inhibiting.csv'))

    def draw_thermodynamics_cdf(self, linewidth=2):
        """
            Test the hypothesis that irreversible reactions are more likely
            to be regulated allosterically
        """
        def comparative_cdf(x, y, data, ax):
            xvals = data[x].unique()
            for xval in xvals:
                d = data.loc[data[x] == xval, y]
                sns.kdeplot(d,
                            cumulative=True, ax=ax, bw=.15,
                            linewidth=linewidth,
                            label=xval + ' (N = %d)' % d.shape[0])
            ax.set_ylim(0, 1)
            ax.set_ylabel(r'Cumulative distribution')
            ax.set_xlabel(y)
            ax.plot([1e3, 1e3], [0, 1], '-', alpha=0.3)

            if len(xvals) == 2:
                ranksum_res = ranksums(data.loc[data[x] == xvals[0], y],
                                       data.loc[data[x] == xvals[1], y])
                ax.set_title(x + '\n$p_{ranksum}$ < %.1g' % ranksum_res.pvalue)
            else:
                ax.set_title(x)

        # get the irreversibility constant (absolute) for all reactions
        # in the BiGG iJO1336 model
        thermo_df = pd.DataFrame.from_csv(settings.ECOLI_THERMO_CACHE_FNAME)

        # remove data about reactions with std=0 (i.e. known values)
        # and reactions with std > 20 (high uncertainty)
        thermo_df = thermo_df[(thermo_df["dG0_prime_std"] > 0) & (thermo_df["dG0_prime_std"] < 20)]

        # select the median value of log(gamma) for each EC number
        # (in general, there should be only one value for each EC number anyway)
        irr_index_l = r"$| log(\Gamma) |$"
        thermo_df[irr_index_l] = thermo_df['logRI'].abs()
        thermo_df = thermo_df[~pd.isnull(thermo_df.EC_number)]

        reg_thermo_df = thermo_df.reset_index()[['EC_number', 'subsystem']].drop_duplicates()
        reg_thermo_df = reg_thermo_df.join(thermo_df.groupby('EC_number').median()[irr_index_l], on='EC_number')

        # print the regulation table with the new irreversibility values
        #reg_thermo_df = self.regulation.join(thermo_df[irr_index_l], on='EC_number')
        #reg_thermo_df.to_csv(os.path.join(settings.RESULT_DIR, 'regulation_with_thermo.csv'))

        # for each EC number, check if it regulated, and if it it positive (+),
        # negative (-) or both (+/-)
        mode_df = self.regulation.groupby('EC_number')['Mode'].apply(set).str.join('/')
        ecs_with_ki = self.regulation.loc[~pd.isnull(self.regulation['KI_Value']), 'EC_number'].unique()

        reg_thermo_df = reg_thermo_df.join(mode_df, on='EC_number')

        reg_thermo_df['Regulation'] = 'regulated'
        reg_thermo_df.loc[pd.isnull(reg_thermo_df['Mode']), 'Regulation'] = 'not regulated'

        reg_thermo_df['Activation'] = 'not activated'
        reg_thermo_df.loc[reg_thermo_df['Mode'].isin(['+', '+/-']), 'Activation'] = 'activated'

        reg_thermo_df['Inhibition'] = 'not inhibited'
        reg_thermo_df.loc[reg_thermo_df['Mode'].isin(['-', '+/-']), 'Inhibition'] = 'inhibited'

        reg_thermo_df['KI_Value'] = 'no $K_I$'
        reg_thermo_df.loc[reg_thermo_df['EC_number'].isin(ecs_with_ki), 'KI_Value'] = 'has $K_I$'

        fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True)
        for i, x_l in enumerate(reg_thermo_df.columns[4:]):
            comparative_cdf(x=x_l, y=irr_index_l,
                            data=reg_thermo_df, ax=axs[i])
            axs[i].set_xlim(0, 15)

        fig.tight_layout()
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_histogram.svg'))
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_histogram.png'), dpi=600)

        # repeat analysis only for CCM reactions
        ccm_thermo_df = reg_thermo_df[
            reg_thermo_df.subsystem.isin(settings.CCM_SUBSYSTEMS)]
        fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True)

        for i, x_l in enumerate(reg_thermo_df.columns[4:]):
            comparative_cdf(x=x_l, y=irr_index_l,
                            data=ccm_thermo_df, ax=axs[i])
            axs[i].set_xlim(0, 15)

        fig.tight_layout()
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_histogram_ccm.svg'))
        fig.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_histogram_ccm.png'), dpi=600)

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

        fig2, axs2 = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
        thermo_df_nz = thermo_df[thermo_df['Num_Regs'] != 0].copy()
        thermo_df_nz['# References / # Regulators'] = \
            thermo_df_nz['Num_Refs'] / thermo_df_nz['Num_Regs']
        thermo_df_nz.plot(kind='scatter', y=irr_index_l,
                          x='# References / # Regulators', ax=axs2[0])
        sns.boxplot(x='Num_Refs', y=irr_index_l, data=thermo_df, ax=axs2[1])
        sns.boxplot(x='Num_Regs', y=irr_index_l, data=thermo_df, ax=axs2[2])
        fig2.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_literature_plot.svg'))
        fig2.savefig(os.path.join(
            settings.RESULT_DIR, 'gibbs_literature_plot.png'), dpi=600)

        return thermo_df

###############################################################################
if __name__ == "__main__":

    #fp = FigurePlotter(rebuild_cache=True)
    fp = FigurePlotter()
    fp.draw_thermodynamics_cdf()

    fp.draw_pathway_met_histogram()
    fp.draw_pathway_histogram()
    fp.draw_venn_diagrams()

    fp.draw_cdf_plots()
    fp.draw_2D_histograms()

    fp.draw_agg_heatmaps(agg_type='gmean')
    fp.draw_agg_heatmaps(agg_type='median')

    fp.draw_full_heapmats()
    fp.draw_full_heapmats(filter_using_model=False)

    fp.print_ccm_table()
    
    plt.close('all')
