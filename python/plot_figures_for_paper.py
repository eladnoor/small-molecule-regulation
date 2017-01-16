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
import os
import json
import seaborn as sns
import numpy as np
from scipy.stats import gmean, ranksums
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import matplotlib
import pdb  # this is a reminder for Elad not to remove this pdb import
from topology import calculate_distances

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

N_ACT_LABEL = 'Number of activating interactions'
N_INH_LABEL = 'Number of inhibiting interactions'

CONDITIONS = ['Glucose', 'Fructose', 'Galactose', 'Gluconate', 'Mannitol',
              'Sorbitol',  'Mannose', 'Glycerol', 'Pyruvate', 'Lactate',
              'Acetate', 'Succinate', 'glcNAc']

class FigurePlotter(object):

    def __init__(self, rebuild_cache=False):
        self.stat_df = pd.DataFrame(index=STAT_TABLE_INDEX,
                                    columns=['km', 'KM_Value',
                                             'regulation', 'KI_Value'])
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

        k['saturation'] = k['concentration'] / (k['concentration'] +
                                                k[value_col])
        k['met:EC'] = k['bigg.metabolite'].str.cat(k['EC_number'], sep=':')
        return k

    @staticmethod
    def calc_agg_sat(k, agg_type='gmean'):
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

        fc_med = fc_med[['saturation']].reset_index()
        fc_med = fc_med.pivot('bigg.metabolite', 'growth condition',
                              'saturation')
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
        self.km = self.km.join(ec2bigg, on='EC_number', how='left')

        self.ki = FigurePlotter.calc_sat(
            self.regulation[~pd.isnull(self.regulation['KI_Value'])],
            'KI_Value', self.met_conc_mean)

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

    def draw_agg_heatmaps(self, agg_type='median'):
        """
            draw heat maps of the [S]/Ki and [S]/Km values across
            the 8 conditions
        """
        def count_values(k, value_col):
            _tmp = k.groupby(('bigg.metabolite', 'growth condition'))
            _tmp = _tmp.count().reset_index().groupby('bigg.metabolite').max()
            return _tmp[value_col].apply(int)

        km_sat_agg = FigurePlotter.calc_agg_sat(self.km, agg_type)
        ki_sat_agg = FigurePlotter.calc_agg_sat(self.ki, agg_type)

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

        fig, ax = plt.subplots(1, 1, figsize=(18, 10))
        clb = matplotlib.colorbar.make_axes(ax)

        sns.heatmap(sat_joined,
                    ax=ax, mask=sat_joined.isnull(), annot=True,
                    cbar=True, vmin=0, vmax=1, cmap='viridis', cbar_ax=clb[0],
                    annot_kws={'fontdict': {'fontsize': 12}})

        # change xtick labels back to the original strings
        # (without the suffixes) and increase the font size
        ax.set_xticklabels(list(km_sat_agg.columns) + list(ki_sat_agg.columns),
                           rotation=90, fontsize=12)

        # rotate the metabolite names back to horizontal, and increase
        # the font size
        ax.set_yticklabels(reversed(sat_joined.index), rotation=0, fontsize=12)

        ax.set_xlabel('growth condition', fontsize=16)
        ax.set_ylabel('')
        ax.set_title('substrates' + ' '*30 + 'inhibitors', fontsize=20)
        clb[0].set_ylabel('%s saturation over all reactions' % agg_type,
                          fontsize=16)
        clb[0].set_yticklabels(np.linspace(0.0, 1.0, 6), fontsize=12)

        ax.axvline(sat_joined.shape[1]/2, 0, 1, color='r')

        settings.savefig(fig, 'heatmap_saturation_%s' % agg_type, dpi=300)

    def draw_full_heapmats(self, filter_using_model=True):
        def pivot_and_sort(k):
            k_piv = k.pivot('met:EC', 'growth condition', 'saturation')
            ind = k_piv.mean(axis=1).sort_values(axis=0, ascending=True).index
            k_piv = k_piv.reindex_axis(ind, axis=0)
            return k_piv

        km = self.km if filter_using_model else self.km_unfiltered
        ki = self.ki if filter_using_model else self.ki_unfiltered
        km_pivoted = pivot_and_sort(km)
        ki_pivoted = pivot_and_sort(ki)
        km_pivoted.index = km_pivoted.index.str.upper()
        ki_pivoted.index = ki_pivoted.index.str.upper()

        # keep and reorder only the conditions that were pre-selected
        km_pivoted = km_pivoted.loc[:, CONDITIONS]
        ki_pivoted = ki_pivoted.loc[:, CONDITIONS]

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
                    cbar=True, vmin=0, vmax=1, annot=True, cmap='viridis',
                    cbar_ax=clb1[0])
        ax1.set_xticklabels(list(ki_pivoted.columns), fontsize=12, rotation=90)
        ax1.set_title('inhibitors', fontsize=20)
        ax1.set_yticklabels(reversed(ki_pivoted.index),
                            rotation=0, fontsize=10)
        ax1.set_xlabel('growth condition', fontsize=16)
        ax1.set_ylabel('')
        clb1[0].set_ylabel('saturation', fontsize=16)
        clb1[0].set_yticklabels(np.linspace(0.0, 1.0, 6), fontsize=12)

        if filter_using_model:
            settings.savefig(fig, 'heatmap_saturation', dpi=200)
            km_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
                                           'heatmap_km_saturation.csv'))
            ki_pivoted.to_csv(os.path.join(settings.RESULT_DIR,
                                           'heatmap_ki_saturation.csv'))
        else:
            settings.savefig(fig, 'heatmap_saturation_unfiltered', dpi=100)

    def draw_venn_diagrams(self):

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
        settings.savefig(fig, 'venn')

        res = {'inhibitors': list(inh_met), 'activators': list(act_met),
               'all_metabolites': list(self.native_mets),
               'inhibited': list(inh_ec), 'activated': list(act_ec),
               'all_reactions': list(self.native_ECs)}
        _fname = os.path.join(settings.RESULT_DIR, 'venn_groups.json')
        with open(_fname, 'w') as fp:
            json.dump(res, fp, indent=4)
        _fname = os.path.join(settings.RESULT_DIR, 'ecoli_interactions.csv')
        self.regulation.to_csv(_fname, 'w')

    def draw_cdf_plots(self, linewidth=2):
        """
            Compare the CDFs of the two fold-change types (for Ki and Km)
        """

        met_color = '#fedf08'  # dandelion
        ki_color = '#fe86a4'  # rosa
        km_color = '#3eaf76'  # dark seafoam green

        fig, axs = plt.subplots(1, 3, figsize=(7.5, 3), sharey=True)

        met_intersection = set(self.km['bigg.metabolite']).intersection(
            self.ki['bigg.metabolite'])
        km_inter = self.km[self.km['bigg.metabolite'].isin(met_intersection)]
        ki_inter = self.ki[self.ki['bigg.metabolite'].isin(met_intersection)]

        ax = axs[0]

        concentrations = pd.melt(self.met_conc_mean)['value']
        concentrations = concentrations[~pd.isnull(concentrations)]

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
        km_values = km_inter.groupby(('met:EC')).first()['KM_Value']
        ki_values = ki_inter.groupby(('met:EC')).first()['KI_Value']
        sns.kdeplot(np.log10(km_values), cumulative=True,
                    ax=ax, bw=.15, color=km_color,
                    label='substrates (N = %d)' % km_values.shape[0],
                    linewidth=linewidth)
        sns.kdeplot(np.log10(ki_values), cumulative=True,
                    ax=ax, bw=.15, color=ki_color,
                    label='inhibitors (N = %d)' % ki_values.shape[0],
                    linewidth=linewidth)
        ax.set_xlim(-2.1, 2.7)
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
        ax.set_ylim(0, 1)
        ax.set_xlabel(r'$K_S$ (in mM)')
        ax.set_title(r'Measured $K_{\rm S}$ values')

        ranksum_res = ranksums(km_values, ki_values)
        ax.text(0.5, 0.3, '$p_{ranksum}$ < %.1g' % ranksum_res.pvalue,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)
        ax.legend(loc='lower right')

        # compare Km and Ki for the intersection of EC numbers

        ax = axs[2]
        ki_saturation = ki_inter['saturation']
        ki_saturation = ki_saturation[~pd.isnull(ki_saturation)]
        km_saturation = km_inter['saturation']
        km_saturation = km_saturation[~pd.isnull(km_saturation)]
        sns.kdeplot(km_saturation, cumulative=True, ax=ax, bw=.01,
                    label='substrates (N = %d)' % km_saturation.shape[0],
                    linewidth=linewidth, color=km_color)
        sns.kdeplot(ki_saturation, cumulative=True, ax=ax, bw=.01,
                    label='inhibitors(N = %d)' % ki_saturation.shape[0],
                    linewidth=linewidth, color=ki_color)

        # Find positions for horizontal annotations
        highval = 0.8
        lowval = 0
        ki_high = float(ki_saturation[ki_saturation <= highval ].shape[0])/ ki_saturation.shape[0]
        ki_low = float(ki_saturation[ki_saturation <= lowval].shape[0])/ ki_saturation.shape[0]

        km_high = float(km_saturation[km_saturation <= highval].shape[0])/ km_saturation.shape[0]
        km_low = float(km_saturation[km_saturation <= lowval].shape[0])/ km_saturation.shape[0]

        # Add vertical lines
        #ax.plot( (lowval,lowval),(0,np.max([ki_low,km_low])),'k--' )
        ax.plot( (highval,highval),(0,np.max([ki_high,km_high])),'k--' )

        # Add horizontal lines
        ax.plot( (1,lowval),(km_low,km_low),color = km_color,linestyle = '--' )
        ax.plot( (1,highval),(km_high,km_high),color = km_color,linestyle = '--' )
        ax.plot( (0,lowval),(ki_low,ki_low),color = ki_color,linestyle = '--' )
        ax.plot( (0,highval),(ki_high,ki_high),color = ki_color,linestyle = '--' )

        # Annotate
        ax.annotate(s='', xytext=(.05,ki_low), xy=(.05,ki_high), arrowprops=dict(facecolor=ki_color, width = 3))

        ax.annotate(s='', xytext=(.85,km_low), xy=(.85,km_high), arrowprops=dict(facecolor=km_color, width = 3))

        ax.text(0.88, 0.2, format((km_high-km_low)*100,'.0f') + '%', horizontalalignment='left', verticalalignment='top',transform=ax.transAxes, color = km_color)

        ax.text(0.1, 0.4, format((ki_high-ki_low)*100,'.0f') + '%', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,color = ki_color)


        ax.grid(visible=False)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0, 1)
        ax.set_xlabel(SAT_FORMULA_S)
        ax.set_title(r'Saturation levels')
        ax.legend(loc='upper left')

        ranksum_res = ranksums(km_saturation, ki_saturation)
        ax.text(0.05, 0.8, '$p_{ranksum}$ < 10$^{%d}$' %
                np.ceil(np.log10(ranksum_res.pvalue)),
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)
        fig.tight_layout()

        settings.savefig(fig, 'cdf_saturation')

    def get_grouped_data(self):
        """
            generate dataframes for metabolite and reactions
            containing counts of activations and inhibitions
        """
        # join interaction table with bigg.reaction IDs
        # and keep only one copy of each reaction-metabolite pair
        cols = ('bigg.reaction', 'bigg.metabolite')
        bigg_effectors = self.regulation.groupby(cols).first().reset_index()
        bigg_effectors.drop('KI_Value', axis=1, inplace=True)

        # add columns for counting the positive and negative interactions
        bigg_effectors[N_ACT_LABEL] = 0
        bigg_effectors[N_INH_LABEL] = 0
        bigg_effectors.loc[bigg_effectors['Mode'] == '+', N_ACT_LABEL] = 1
        bigg_effectors.loc[bigg_effectors['Mode'] == '-', N_INH_LABEL] = 1

        grouped_by_met = bigg_effectors.groupby('bigg.metabolite').sum()
        grouped_by_rxn = bigg_effectors.groupby('bigg.reaction').sum()
        return grouped_by_met, grouped_by_rxn

    def draw_2D_histograms(self, tax2use='kingdom', minsize=10):
        """
            Draw 2D histograms of the number of activating and
            inhibiting reactions
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
            sns.heatmap(np.log2(H.T), mask=(H.T == 0), ax=ax, annot=False,
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
                    # print i, x_i, y_i, H[x_i, y_i]
                    ax.annotate(i, xy=(x_i, ymax-y_i),
                                xytext=(x_i+1, ymax-y_i-1),
                                ha='center', va='top', size=6,
                                textcoords='data')

        def plot_jointhist(data, xlabel, ylabel, xmax, ymax):
            """
                plot the histogram as a scatter plot with marginal histograms,
                and ensure empty bins are easily distinguishable from ones
                that have at least 1 hit.
            """
            x = data[xlabel]
            y = data[ylabel]

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
            ann_df = ann_df[(ann_df[xlabel] > 12) | (ann_df[ylabel] > 5)]
            for i, row in ann_df.iterrows():
                    plt.annotate(i, xy=(row[xlabel], row[ylabel]),
                                 xytext=(row[xlabel]+1, row[ylabel]+1),
                                 ha='center', va='top', size=10,
                                 textcoords='data')

            # Next, plot the marginal histograms
            g = g.plot_marginals(sns.distplot, kde=False)

            return plt.gcf()

        grouped_by_met, grouped_by_rxn = self.get_grouped_data()

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
        xmax = max(grouped_by_met[N_INH_LABEL].max(),
                   grouped_by_rxn[N_INH_LABEL].max())
        ymax = max(grouped_by_met[N_ACT_LABEL].max(),
                   grouped_by_rxn[N_ACT_LABEL].max())

        plot_2d_hist(grouped_by_met, N_INH_LABEL, N_ACT_LABEL,
                     axs[0], xmax, ymax)
        axs[0].set_xlabel('')
        plot_2d_hist(grouped_by_rxn, N_INH_LABEL, N_ACT_LABEL,
                     axs[1], xmax, ymax)

        settings.savefig(fig, '2D_histograms_heatmap')

        fig = plot_jointhist(grouped_by_met, N_INH_LABEL, N_ACT_LABEL,
                             xmax, ymax)
        settings.savefig(fig, 'jointplot_met')

        fig = plot_jointhist(grouped_by_rxn, N_INH_LABEL, N_ACT_LABEL,
                             xmax, ymax)
        settings.savefig(fig, 'jointplot_rxn')

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

    def draw_pathway_histogram(self):
        pathways = set(self.regulation['bigg.subsystem.reaction'])
        pathways.update(self.regulation['bigg.subsystem.metabolite'])
        pathways = sorted(pathways)
        if np.nan in pathways:
            pathways.remove(np.nan)

        cols = ['bigg.metabolite', 'bigg.subsystem.metabolite',
                'bigg.reaction', 'bigg.subsystem.reaction', 'Mode']
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
                        hist_mat_act[i, j] = act_table.at[pathways[i],
                                                          pathways[j]]

        N = len(pathways)
        hist_mat_inh = np.zeros((N, N))
        for i in xrange(N):
            if pathways[i] in inh_table.index:
                for j in xrange(N):
                    if pathways[j] in inh_table.columns:
                        hist_mat_inh[i, j] = inh_table.at[pathways[i],
                                                          pathways[j]]

        vmax = max(hist_mat_act.max(), hist_mat_inh.max())

        fig, axs = plt.subplots(1, 2, figsize=(20, 9))

        ax = axs[0]
        sns.heatmap(hist_mat_act, ax=ax, annot=False, cbar=False,
                    cmap='viridis', xticklabels=pathways, yticklabels=pathways,
                    vmin=0, vmax=vmax)
        ax.set_title('activation')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('activated enzyme subsystem')

        ax = axs[1]
        sns.heatmap(hist_mat_inh, ax=ax, annot=False, cbar=True,
                    cmap='viridis', xticklabels=pathways, yticklabels=False,
                    vmin=0, vmax=vmax)
        ax.set_title('inhibition')
        ax.set_xlabel('inhibited enzyme subsystem')

        fig.tight_layout(pad=4)
        settings.savefig(fig, 'pathway_histograms', dpi=300)
        act_table.to_csv(os.path.join(settings.RESULT_DIR,
                                      'pathway_histograms_activating.csv'))
        inh_table.to_csv(os.path.join(settings.RESULT_DIR,
                                      'pathway_histograms_inhibiting.csv'))

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

        act_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_activating.csv'))
        inh_table.to_csv(os.path.join(
            settings.RESULT_DIR, 'pathway_met_histograms_inhibiting.csv'))

        # Remove metabolites and pathways with insufficient numbers of
        # data points
        act2plot = act_table[act_table.sum(axis=1).fillna(0) > 3]
        dropcol = act2plot.columns[act2plot.sum(axis=0).fillna(0) < 3]
        act2plot = act2plot.drop(dropcol, axis=1)

        inh2plot = inh_table[inh_table.sum(axis=1).fillna(0) > 7]
        dropcol = inh2plot.columns[inh2plot.sum(axis=0).fillna(0) < 7]
        inh2plot = inh2plot.drop(dropcol, axis=1)

        # Plot activating and inhibiting matrices in seaborn

        act_roworder = FigurePlotter.cluster_matrix(act2plot)
        act_colorder = FigurePlotter.cluster_matrix(act2plot.T)

        inh_roworder = FigurePlotter.cluster_matrix(inh2plot)
        inh_colorder = FigurePlotter.cluster_matrix(inh2plot.T)

        fig, axs = plt.subplots(1, 2, figsize=(10, 10), sharey=False)

        ax = axs[0]
        sns.heatmap(act2plot.ix[act_roworder, act_colorder],
                    ax=ax, annot=False, cbar=True, cmap='viridis')
        ax.set_title('activation')
        ax.set_ylabel(ylabel)
        ax.set_xlabel('activated enzyme subsystem')
        for label in ax.get_xticklabels():
            label.set(rotation=90)
        for label in ax.get_yticklabels():
            label.set(rotation=0)

        ax = axs[1]
        sns.heatmap(inh2plot.ix[inh_roworder, inh_colorder],
                    ax=ax, annot=False, cbar=True, cmap='viridis')
        ax.set_title('inhibition')
        ax.set_xlabel('inhibited enzyme subsystem')
        for label in ax.get_xticklabels():
            label.set(rotation=90)
        for label in ax.get_yticklabels():
            label.set(rotation=0)

        fig.tight_layout(pad=4)
        settings.savefig(fig, 'pathway_met_histograms', dpi=300)

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
    def comparative_cdf(x, y, data, ax, linewidth=2):
        xvals = data[x].unique()
        for xval in xvals:
            d = data.loc[data[x] == xval, y]
            sns.kdeplot(d, cumulative=True, ax=ax, bw=.15,
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

    def draw_thermodynamics_cdf(self):
        """
            Test the hypothesis that irreversible reactions are more likely
            to be regulated allosterically
        """
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
        irr_index_l = r"$| log(\Gamma) |$"
        thermo_df[irr_index_l] = thermo_df['logRI'].abs()
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

        # for each EC number, check if it regulated, and if it it positive (+),
        # negative (-) or both (+/-)
        mode_df = self.regulation.groupby('EC_number')['Mode']
        mode_df = mode_df.apply(set).str.join('/')  # convert to short string
        ecs_with_ki = self.regulation.loc[
            ~pd.isnull(self.regulation['KI_Value']), 'EC_number'].unique()

        reg_thermo_df = reg_thermo_df.join(mode_df, on='EC_number')

        reg_thermo_df['Regulation'] = 'regulated'
        unreg_ind = pd.isnull(reg_thermo_df['Mode'])
        reg_thermo_df.loc[unreg_ind, 'Regulation'] = 'not regulated'

        reg_thermo_df['Activation'] = 'not activated'
        act_ind = reg_thermo_df['Mode'].isin(['+', '+/-'])
        reg_thermo_df.loc[act_ind, 'Activation'] = 'activated'

        reg_thermo_df['Inhibition'] = 'not inhibited'
        inh_ind = reg_thermo_df['Mode'].isin(['-', '+/-'])
        reg_thermo_df.loc[inh_ind, 'Inhibition'] = 'inhibited'

        reg_thermo_df['KI_Value'] = 'no $K_I$'
        ki_ind = reg_thermo_df['EC_number'].isin(ecs_with_ki)
        reg_thermo_df.loc[ki_ind, 'KI_Value'] = 'has $K_I$'

        fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True)
        for i, x_l in enumerate(reg_thermo_df.columns[4:]):
            FigurePlotter.comparative_cdf(x=x_l, y=irr_index_l,
                                          data=reg_thermo_df, ax=axs[i])
            axs[i].set_xlim(0, 15)

        fig.tight_layout()
        settings.savefig(fig, 'gibbs_histogram')

        # repeat analysis only for CCM reactions
        ccm_thermo_df = reg_thermo_df[
            reg_thermo_df.subsystem.isin(settings.CCM_SUBSYSTEMS)]
        fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True)

        for i, x_l in enumerate(reg_thermo_df.columns[4:]):
            FigurePlotter.comparative_cdf(x=x_l, y=irr_index_l,
                                          data=ccm_thermo_df, ax=axs[i])
            axs[i].set_xlim(0, 15)

        fig.tight_layout()
        settings.savefig(fig, 'gibbs_histogram_ccm')

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

        fig2, axs2 = plt.subplots(1, 4, figsize=(14, 4), sharey=True)
        thermo_df_nz = thermo_df[thermo_df['Num_Regs'] != 0].copy()
        thermo_df_nz['# References / # Regulators'] = \
            thermo_df_nz['Num_Refs'] / thermo_df_nz['Num_Regs']
        thermo_df_nz.plot(kind='scatter', y=irr_index_l,
                          x='# References / # Regulators', ax=axs2[0])
        sns.boxplot(x='Num_Refs', y=irr_index_l, data=thermo_df, ax=axs2[1])
        sns.boxplot(x='Num_Regs', y=irr_index_l, data=thermo_df, ax=axs2[2])

        from scipy.stats import mannwhitneyu as mwu
        regulated = thermo_df.ix[thermo_df['is regulated'] == 'Yes',
                                 irr_index_l]
        notregulated = thermo_df.ix[thermo_df['is regulated'] == 'No',
                                    irr_index_l]
        stat, p = mwu(regulated, notregulated)

        sns.boxplot(x='is regulated', y=irr_index_l, data=thermo_df,
                    ax=axs2[3])
        axs2[3].set_title('Mann Whitney P Value = ' + str(p))
        settings.savefig(fig2, 'gibbs_literature_plot')

        return thermo_df

    def draw_ccm_thermodynamics_cdf(self):
        ccm_thermo_df = pd.DataFrame.from_csv(settings.ECOLI_CCM_THERMO_FNAME)
        irr_index_l = r"$| log(\Gamma) |$"
        ccm_thermo_df[irr_index_l] = ccm_thermo_df['logRI'].abs()

        from scipy.stats import mannwhitneyu as mwu
        fig, ax = plt.subplots(1, 1, figsize=(8, 8), sharey=True)
        sns.boxplot(x='is regulated', y=irr_index_l, data=ccm_thermo_df, ax=ax)
        plt.xlabel('Is Regulated')
        plt.ylabel('Log Reversibility Index')

        regulated = ccm_thermo_df.ix[ccm_thermo_df['is regulated'] == 'yes',
                                     irr_index_l]
        notregulated = ccm_thermo_df.ix[ccm_thermo_df['is regulated'] == 'no',
                                        irr_index_l]
        stat, p = mwu(regulated, notregulated)

        plt.title('Mann Whitney P Value: ' + str(p))
        settings.savefig(fig, 'gibbs_boxplot_ccm_curated')
        plt.close('boxplot')

        fig, ax = plt.subplots(1, 1, figsize=(3, 3))
        FigurePlotter.comparative_cdf(x='is regulated', y=irr_index_l,
                                      data=ccm_thermo_df, ax=ax)
        ax.set_xlabel(irr_index_l)
        ax.set_xlim(0, 15)

        fig.tight_layout()
        settings.savefig(fig, 'gibbs_histogram_ccm_curated')

    def compare_km_ki(self, filter_using_model=False):

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

        res['KI_Values'] = ki.groupby('bigg.metabolite')['KI_Value'].mean()
        res['KI_Number'] = ki.groupby('bigg.metabolite')['EC_number'].nunique()
        res['KM_Values'] = km.groupby('bigg.metabolite')['KM_Value'].mean()
        res['KM_Number'] = km.groupby('bigg.metabolite')['EC_number'].nunique()

        # Drop rows where we don't have data for both
        res = res.dropna()

        # Keep only metabolites with at least 2 measurements of each
        res = res[res['KI_Number'] > 1]
        res = res[res['KM_Number'] > 1]

        res['PValue'] = np.nan

        # for each metabolite, if there is sufficient data, test
        for ii in res.index:
            kid = ki[ki['bigg.metabolite'] == ii]['KI_Value']
            kmd = km[km['bigg.metabolite'] == ii]['KM_Value']

            s,p = st.mannwhitneyu( kid,kmd )
            res.at[ii,'PValue'] = p
            res['QValue'] = padjust(res['PValue'],method = 'fdr_bh')[1]
        res = res.sort_values('PValue')

        maxval = 2*np.max( [res['KI_Values'].max(),res['KM_Values'].max()] )
        minval = 0.5*np.min( [res['KI_Values'].min(),res['KM_Values'].min()] )

        fig,ax = plt.subplots(figsize = (8,8))
        ax.scatter(res['KI_Values'],res['KM_Values'],s = 10*res['KI_Number'], color = 'grey')
        ax.axis([minval,maxval,minval,maxval])
        ax.set_xscale('log')
        ax.set_yscale('log')

        # Calculate log ratio of values
        res['Ratio'] = np.log2( res['KI_Values'] / res['KM_Values'] )

        for ii in res.index:
            if res.at[ii,'QValue'] < 0.1:
            	ax.scatter(res.at[ii,'KI_Values'], res.at[ii,'KM_Values'], color = 'r', s = 10*res.at[ii,'KI_Number'])
                ax.text(1.2*res.at[ii,'KI_Values'],res.at[ii,'KM_Values'],ii)
        plt.xlabel('Mean KI')
        plt.ylabel('Mean KM')

        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

        settings.savefig(fig, 'km_vs_ki')

        # Make volcano plot
        fig,ax = plt.subplots(figsize = (8,8))
        ax.scatter(res['Ratio'],-np.log10(res['QValue']),color = 'grey',s = 4*(res['KI_Number'] + res['KM_Number']) )

        for ii in res.index:
            if res.at[ii,'QValue'] < 0.1 and np.abs(res.at[ii,'Ratio']) > 1:
            	ax.scatter(res.at[ii,'Ratio'], -np.log10(res.at[ii,'QValue']), color = 'r', s = 4*(res['KI_Number'] + res['KM_Number']) )
                ax.text(1.2*res.at[ii,'Ratio'], -np.log10(res.at[ii,'QValue']),ii)

        plt.xlabel('Log2 (Mean KI/Mean KM)')
        plt.ylabel('-Log10 Q Value')

        # Plot some lines
        highq = -np.log10(res['QValue'].min()) + .1
        highfc = res['Ratio'].max() + 1
        ax.axis([-highfc,highfc,0,highq])
        #ax.plot( (0,0),(0,highq),'k--' ) # vertical line
        #ax.plot( (-highfc,highfc),(0,0),'k--' ) # horizontal line
        settings.savefig(fig, 'km_vs_ki_volcano')

    def draw_distance_histograms(self):
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
        axs[0, 0].hist(all_distances, color='#939598', **args)
        axs[1, 0].hist(smrn_dist['distance'], color='#8a5ec9', **args)
        axs[0, 1].hist(inh_distances, color='#f5821f', **args)
        axs[1, 1].hist(act_distances, color='#00aeef', **args)

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

        fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.pdf'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_distances.png'))
        smrn_dist.to_csv(os.path.join(settings.CACHE_DIR,
                                      'iJO1366_SMRN_dist.csv'), index=False)

    def draw_degree_histograms(self):
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
        axs[0, 0].hist(rxn_hist['bigg.metabolite'], color='#808080', **args)
        axs[1, 0].hist(met_hist['bigg.reaction'], color='#808080', **args)
        axs[0, 1].hist(rxn_hist_ccm['bigg.metabolite'], color='#7060ef', **args)
        axs[1, 1].hist(met_hist_ccm['bigg.reaction'], color='#7060ef', **args)

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
        fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_degrees.pdf'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'SMRN_degrees.png'),
                    dpi=300)

###############################################################################
if __name__ == "__main__":
    plt.close('all')
#    fp = FigurePlotter(rebuild_cache=True)
    fp = FigurePlotter()
    fp.draw_2D_histograms()
    fp.draw_thermodynamics_cdf()

    fp.draw_ccm_thermodynamics_cdf()

    fp.draw_pathway_met_histogram()
    fp.draw_pathway_histogram()
    fp.draw_venn_diagrams()

    fp.draw_cdf_plots()

    fp.draw_agg_heatmaps(agg_type='gmean')
    fp.draw_agg_heatmaps(agg_type='median')

    fp.draw_full_heapmats()
    fp.draw_full_heapmats(filter_using_model=False)

    fp.print_ccm_table()
    fp.compare_km_ki()

    fp.draw_degree_histograms()
#    fp.draw_distance_histograms()

    plt.close('all')
