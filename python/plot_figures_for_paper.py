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
import os, re, json
import seaborn as sns
import numpy as np
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import matplotlib
sns.set('paper', style='white')

ORGANISM = 'Escherichia coli'
SAT_FORMULA_S = r'$[S]/\left([S] + K_S\right)$'
SAT_FORMULA_M = r'$[S]/\left([S] + K_M\right)$'
SAT_FORMULA_I = r'$[S]/\left([S] + K_I\right)$'

class FigurePlotter(object):
    
    def __init__(self, rebuild_cache=False):
        self.bigg = BiGG()
        self.kegg = KEGG()
        
        # load the BiGG model for E. coli in order to define the full scope
        # of metabolites and reactions in this organism
        self.model, self.S = settings.get_ecoli_json()    
        
        if rebuild_cache:
            map_ligands.rebuild_cache()
        
        self.get_data()
    
    def draw_sankey_diagram(self):
        all_brenda_interactions = 100000 # both activation and inhibition
        
        # TODO: this is just a fake figure. Needs to be replaced with real data
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sankey = Sankey(ax=ax, scale=0.01, offset=0.2, head_angle=180,
                        format='%.0f', unit='%')
        sankey.add(flows=[25, 0, 60, -10, -20, -5, -15, -10, -40],
                   labels=['', '', '', 'First', 'Second', 'Third', 'Fourth',
                           'Fifth', 'Hurray!'],
                   orientations=[-1, 1, 0, 1, 1, 1, -1, -1, 0],
                   pathlengths=[0.25, 0.25, 0.25, 0.25, 0.25, 0.6, 0.25, 0.25,
                                0.25],
                   patchlabel="Widget\nA",
                   alpha=0.2, lw=2.0)  # Arguments to matplotlib.patches.PathPatch()
        diagrams = sankey.finish()
        diagrams[0].patch.set_facecolor('#37c959')
        diagrams[0].texts[-1].set_color('r')
        diagrams[0].text.set_fontweight('bold')
        fig.savefig(os.path.join(settings.RESULT_DIR, 'sankey.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'sankey.png'), dpi=600)
    
    @staticmethod
    def get_kinetic_param(name, value_col=None, organism=ORGANISM):
        k = settings.read_cache(name)
        
        # filter by organsim
        k = k[k['Organism'] == organism]
    
        # filter out mutated enzymes
        if 'Commentary' in k.columns:
            k = k[k['Commentary'].str.find('mutant') == -1]
            k = k[k['Commentary'].str.find('mutation') == -1]
        
        # remove values with unmatched ligand
        k = k[pd.notnull(k['bigg.metabolite'])]
        k['bigg.metabolite'] = k['bigg.metabolite'].str.lower()
    
        if value_col is not None:
            # remove entries lacking quantitative data
            k = k[k[value_col] > 0]
            return k[['EC_number', 'bigg.metabolite', value_col]]
        else:
            return k[['EC_number', 'bigg.metabolite']]        
    
    @staticmethod
    def calc_sat(k, value_col, conc_df):
        # choose the minimum value among all repeats    
        k = k.groupby(['EC_number', 'bigg.metabolite'])[value_col].min().reset_index()
        
        # join data with measured concentrations
        k = k.join(conc_df, on='bigg.metabolite', how='inner')
        
        # melt table so each line will be a combination of EC, substrate/inhibitor and
        # growth condition
        k = pd.melt(k, id_vars=('EC_number', 'bigg.metabolite', value_col),
                    var_name='growth condition', value_name='concentration')
        
        k['saturation'] = k['concentration'] / (k['concentration'] + k[value_col])
        k['met:EC'] = k['bigg.metabolite'].str.cat(k['EC_number'], sep=':')
        return k
    
    @staticmethod
    def calc_median_sat(k):
        """
            calculates the [S]/K_S for all matching EC-metabolite pairs,
            in log2-fold-change.
            
            Input:
                K_df    - a DataFrame with three columns: EC_number, bigg.metabolite, Value
                conc_df - a DataFrame with 
        """
        fc_med = k.groupby(('bigg.metabolite', 'growth condition')).median()[['saturation']].reset_index()
        fc_med = fc_med.pivot('bigg.metabolite', 'growth condition', 'saturation')
        return fc_med.sort_index(axis=0)
    
    def get_data(self):
        _df = pd.DataFrame.from_csv(settings.ECOLI_METAB_FNAME)
        _df.index.name = 'bigg.metabolite'
        self.met_conc_mean = _df.iloc[:, 1:9]
        
        # remove the _c suffix from the compound names
        self.met_conc_mean.index = map(lambda s: s[0:-2], self.met_conc_mean.index)
        #met_conc_std = _df.iloc[:, 10:]
        
        colmap = dict(map(lambda x: (x, x[:-7]), self.met_conc_mean.columns))
        self.met_conc_mean.rename(columns=colmap, inplace=True)
        
        km_raw = FigurePlotter.get_kinetic_param('km', 'KM_Value')
        self.km = FigurePlotter.calc_sat(km_raw, 'KM_Value', self.met_conc_mean)
        
        ki_raw = FigurePlotter.get_kinetic_param('ki', 'KI_Value')
        self.ki = FigurePlotter.calc_sat(ki_raw, 'KI_Value', self.met_conc_mean)
    
        self.ki.to_csv(os.path.join(settings.RESULT_DIR, 'ki_vs_conc.csv'))
        self.km.to_csv(os.path.join(settings.RESULT_DIR, 'km_vs_conc.csv'))

        act = FigurePlotter.get_kinetic_param('activating', None)
        inh = FigurePlotter.get_kinetic_param('inhibiting', None)
        act['type'] = 'activation'
        inh['type'] = 'inhibition'
        self.interactions = pd.concat([act, inh])
    
    def draw_median_heatmaps(self):
        """
            draw heat maps of the [S]/Ki and [S]/Km values across the 8 conditions
        """
        km_sat_med = FigurePlotter.calc_median_sat(self.km)
        ki_sat_med = FigurePlotter.calc_median_sat(self.ki)
        
        sat_joined = km_sat_med.join(ki_sat_med, how='inner', lsuffix='_sub', rsuffix='_inh')
        sat_joined = sat_joined.reindex_axis(sat_joined.mean(axis=1).sort_values(axis=0, ascending=False).index, axis=0)
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        clb = matplotlib.colorbar.make_axes(ax)
        
        sns.heatmap(sat_joined, ax=ax, mask=sat_joined.isnull(), annot=True, 
                    cbar=True, vmin=0, vmax=1, cmap='viridis', cbar_ax=clb[0],
                    annot_kws={'fontdict': {'fontsize': 12}})
        
        # change xtick labels back to the original strings (without the suffixes) and increase the font size
        ax.set_xticklabels(list(ki_sat_med.columns) + list(ki_sat_med.columns), rotation=90, fontsize=12)
        
        # rotate the metabolite names back to horizontal, and increase the font size
        ax.set_yticklabels(reversed(sat_joined.index), rotation=0, fontsize=12)
        
        ax.set_xlabel('growth condition', fontsize=16)
        ax.set_ylabel('')
        ax.set_title('substrates' + ' '*30 + 'inhibitors', fontsize=20)
        clb[0].set_ylabel('median saturation over all reactions', fontsize=16)
        clb[0].set_yticklabels(np.linspace(0.0, 1.0, 6), fontsize=12)
        
        ax.axvline(sat_joined.shape[1]/2, 0, 1, color='r')
        
        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation_median.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation_median.png'), dpi=300)
    
    def draw_full_heapmats(self):
    
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 30))
    
        km_pivoted = self.km.pivot('met:EC', 'growth condition', 'saturation')
        km_pivoted = km_pivoted.reindex_axis(km_pivoted.mean(axis=1).sort_values(axis=0, ascending=True).index, axis=0)
        sns.heatmap(km_pivoted, ax=ax0, mask=km_pivoted.isnull(),
                    cbar=False, vmin=0, vmax=1, cmap='viridis')
        ax0.set_xticklabels(list(km_pivoted.columns), fontsize=12, rotation=90)
        ax0.set_yticklabels(reversed(km_pivoted.index), rotation=0, fontsize=6)
        ax0.set_title('substrates', fontsize=20)
        ax0.set_xlabel('growth condition', fontsize=16)
        ax0.set_ylabel('')
    
        ki_pivoted = self.ki.pivot('met:EC', 'growth condition', 'saturation')
        ki_pivoted = ki_pivoted.reindex_axis(ki_pivoted.mean(axis=1).sort_values(axis=0, ascending=True).index, axis=0)
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
        
        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation.svg'))
        fig.savefig(os.path.join(settings.RESULT_DIR, 'heatmap_saturation.png'), dpi=100)
        km_pivoted.to_csv(os.path.join(settings.RESULT_DIR, 'heatmap_km_saturation.csv'))
        ki_pivoted.to_csv(os.path.join(settings.RESULT_DIR, 'heatmap_ki_saturation.csv'))
    
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
        
        # make a list of all the EC numbers which are mapped to a BiGG reaction 
        # in the E. coli model, which in turn is mapped to at least one E. coli gene
        rids_with_genes = set()       
        for d in self.model['reactions']:
            rid = d['id']
            rule = d['gene_reaction_rule']
            if re.search('b[0-9]+', rule) is not None:
                rids_with_genes.add(rid.lower())
        
        # use the self.bigg object to convert these BiGG IDs to EC numbers
        bigg_reactions = self.bigg.reaction_df
        native_ec = set(bigg_reactions[bigg_reactions['bigg.reaction'].isin(rids_with_genes)].index)

        # keep only interactions that are with a native EC number
        native_interactions = self.interactions[self.interactions['EC_number'].isin(native_ec)]

        # make a list of all the metabolites that are cytoplasmic
        mets_in_cytoplasm = set()       
        for d in self.model['metabolites']:
            met = d['id']
            if d['compartment'] == u'c':
                mets_in_cytoplasm.add(met.lower()[:-2])
        
        # keep only interactions that involve a small-molecule that is native
        # in the cytoplasm
        native_interactions = native_interactions[native_interactions['bigg.metabolite'].isin(mets_in_cytoplasm)]
        
        print "found %d native interactions in %s" % (native_interactions.shape[0], ORGANISM)
        
        ind_inh = native_interactions['type']=='inhibition'
        ind_act = native_interactions['type']=='activation'

        inh_met = set(native_interactions.loc[ind_inh, 'bigg.metabolite'])
        act_met = set(native_interactions.loc[ind_act, 'bigg.metabolite'])
        inh_ec = set(native_interactions.loc[ind_inh, 'EC_number'])
        act_ec = set(native_interactions.loc[ind_act, 'EC_number'])
        
        fig, axs = plt.subplots(1, 2, figsize=(7, 5))
        venn3_sets(inh_met, act_met, mets_in_cytoplasm,
                   set_labels=('inhibitors', 'activators', 'E. coli metabolites'), ax=axs[0])
        venn3_sets(inh_ec, act_ec, native_ec,
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
               'all_metabolites': list(mets_in_cytoplasm),
               'inhibited': list(inh_ec), 'activated': list(act_ec),
               'all_reactions': list(native_ec)}
        with open(os.path.join(settings.RESULT_DIR, 'venn_groups.json'), 'w') as fp:
            json.dump(res, fp, indent=4)
        native_interactions.to_csv(open(os.path.join(settings.RESULT_DIR, 'ecoli_interactions.csv'), 'w'))

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
        
        concentrations = pd.melt(fp.met_conc_mean)['value']
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
            

        # join interaction table with bigg.reaction IDs (using EC numbers)
        # and keep only one copy of each reaction-metabolite pair
        bigg_effectors = self.interactions.join(self.bigg.reaction_df,
                                                how='inner', on='EC_number')
        cols = ('bigg.reaction', 'bigg.metabolite')
        bigg_effectors = bigg_effectors.groupby(cols).first().reset_index()
        
        # add columns for counting the positive and negative interactions
        n_act_label = 'Number of activating interactions'        
        n_inh_label = 'Number of inhibiting interactions'
        
        bigg_effectors[n_act_label] = 0
        bigg_effectors[n_inh_label] = 0
        bigg_effectors.loc[bigg_effectors['type'] == 'activation', n_act_label] = 1
        bigg_effectors.loc[bigg_effectors['type'] == 'inhibition', n_inh_label] = 1
        
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

###############################################################################
if __name__ == "__main__":
    
    #fp = FigurePlotter(rebuild_cache=True)
    fp = FigurePlotter()
    
    #fp.draw_sankey_diagram()
    #fp.draw_median_heatmaps()
    #fp.draw_full_heapmats()
    fp.draw_venn_diagrams()
    #fp.draw_cdf_plots()
    fp.draw_2D_histograms()
    
