# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 19:49:04 2016

@author: noore
"""
import os
import pandas as pd
import settings as S
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from bool_parser import BoolParser
import matplotlib.pyplot as plt
import seaborn as sns

class Volume(object):

    def __init__(self):
        self.cobra_model = create_cobra_model_from_sbml_file(S.ECOLI_SBML_FNAME)
        convert_to_irreversible(self.cobra_model)

        self.met_conc_df = self.get_metabolite_data()
        self.km_df = self.get_km_data()
        self.enz_conc_df, self.enz_mw_df = self.get_enzyme_data()
        self.flux_df = self.get_flux_data()
    
        self.data_df = pd.merge(self.km_df, self.met_conc_df, on='bigg.metabolite')
        self.data_df = pd.merge(self.data_df, self.enz_conc_df, on=['bigg.reaction', 'condition'])
        self.data_df = pd.merge(self.data_df, self.enz_mw_df, on=['bigg.reaction'])
        self.data_df = pd.merge(self.data_df, self.flux_df, on=['bigg.reaction', 'condition'])
        
        # keep only rows with non-zero flux, non-zero enzyme, and stoichiometry coeff = 1
        
        ind = (self.data_df['flux [mmol/gCDW/h]'] > 0) & \
              (self.data_df['enzyme conc [M]'] > 0) & \
              (self.data_df['stoichiometry'] == -1)
        
        self.data_df = self.data_df[ind]
        
    def get_metabolite_data(self):
        """
            get the metabolomics data from Gerosa et al. 2015
        """
        _df = pd.DataFrame.from_csv(S.ECOLI_METAB_FNAME)
        _df.index.name = 'bigg.metabolite'
        met_conc_mean = _df.iloc[:, 1:9] * 1e-3 # convert mM to M
        met_conc_std = _df.iloc[:, 10:]  * 1e-3 # convert mM to M
        
        met_conc_mean.columns = [c[:-7].lower() for c in met_conc_mean.columns]
        met_conc_std.columns = [c[:-6].lower() for c in met_conc_std.columns]
        
        met_conc_df = pd.melt(met_conc_mean.reset_index(), id_vars='bigg.metabolite',
                              var_name='condition', value_name='metabolite conc [M]')
    
        # join the concentration data with the molecular weight of each metabolite
        met_mw_df = pd.DataFrame(columns=('bigg.metabolite', 'metabolite MW [Da]'),
            data=[(m.id, m.formula_weight) for m in self.cobra_model.metabolites])
        met_mw_df.set_index('bigg.metabolite', inplace=True)
        met_conc_df = met_conc_df.join(met_mw_df, on='bigg.metabolite')
        
        return met_conc_df
    
    def get_km_data(self):
        km_df = S.read_cache('km')
        km_df = km_df[km_df['Organism'] == 'Escherichia coli']
        km_df = km_df[km_df['KM_Value'] != -999]
        km_df = km_df[['EC_number', 'KM_Value', 'bigg.metabolite']]
        km_df = km_df.groupby(('EC_number', 'bigg.metabolite')).median().reset_index()

        # some compounds have specific steriochemistry in BRENDA, but not in the
        # E. coli model (and other datasets). Therefore, we need to map them to 
        # the stereo-unspecific BiGG IDs in order to join the tables later
        stereo_mapping = {'fdp_B_c': 'fdp_c', 'f6p_B_c': 'f6p_c'}
        km_df['bigg.metabolite'].replace(stereo_mapping, inplace=True)
        
        # get a mapping from EC numbers to bigg.reaction,
        # remember we need to duplicate every reaction ID also for the reverse
        # reaction (since we use a model that is converted to irreversible)
        model_reactions = S.get_reaction_table_from_xls()
        bigg2ec = model_reactions[['Reaction Abbreviation', 'EC Number']]
        bigg2ec.rename(columns={'Reaction Abbreviation': 'bigg.reaction',
                                'EC Number': 'EC_number'}, inplace=True)
        bigg2ec = bigg2ec[~pd.isnull(bigg2ec['EC_number'])]
        bigg2ec['bigg.reaction'] = bigg2ec['bigg.reaction'].str.lower()
       
        bigg2ec_rev = bigg2ec.copy()
        bigg2ec_rev['bigg.reaction'] = bigg2ec_rev['bigg.reaction'].apply(lambda s: s + '_reverse')
        bigg2ec = pd.concat([bigg2ec, bigg2ec_rev], axis=0)
    
        # get the stoichiometric matrix in order to match only the substrates
        # of each reaction (i.e. the products should be associated only with
        # the "reverse" reaction)
        stoich_df = []
        for r in self.cobra_model.reactions:
            for m, coeff in r.metabolites.iteritems():
                stoich_df.append([r.id.lower(), m.id, coeff])
        stoich_df = pd.DataFrame(columns=['bigg.reaction', 'bigg.metabolite', 'stoichiometry'],
                                 data=stoich_df)
    
        km_df = pd.merge(km_df, bigg2ec, on='EC_number')
        km_df = pd.merge(km_df, stoich_df, on=('bigg.reaction', 'bigg.metabolite'))
        km_df['Km [M]'] = km_df['KM_Value'] * 1e-3 # convert mM to M
        km_df.drop('KM_Value', axis=1, inplace=True)
        
        return km_df
    
    def get_enzyme_data(self):
        """
            get the enzyme data from Schmidt et al. 2015
        """
        _df = pd.DataFrame.from_csv(S.ECOLI_PROT_FNAME)
        _df = _df[~pd.isnull(_df['Bnumber'])]
        _df.set_index('Bnumber', inplace=True)
        enz_mw         = _df['Molecular weight (Da)']
        enz_count_mean = _df.iloc[:, 6:28] # counts per 10^-15 liter
        enz_count_cv   = _df.iloc[:, 29:51]
        enz_count_mean.columns = [c[:-7].lower() for c in enz_count_mean.columns]
        enz_count_cv.columns = [c[:-5].lower() for c in enz_count_cv.columns]
        
        enz_count_mean = enz_count_mean.apply(pd.to_numeric).fillna(0.0)
        
        # for each reaction, get a single value for the enzyme concentration
        # in [M] and the molecular weight of the complex
        
        reaction2boolparser = {}
        for r in self.cobra_model.reactions:
            if r.gene_reaction_rule == '':
                continue
            if r.gene_reaction_rule.find('s0001') != -1:
                continue
            reaction2boolparser[r.id.lower()] = BoolParser(r.gene_reaction_rule)
    
        # calculate the MW of each enzyme complex, using the gene reaction rules
        # from the iJO1336 model    
        enz_mw_df = pd.DataFrame(index=reaction2boolparser.keys(),
                                 columns=['enzyme MW [Da]'])
        enz_mw_df.index.name = 'bigg.reaction'
        enzyme_mw_dict = enz_mw.to_dict()
        for r, bool_parser in reaction2boolparser.iteritems():
            enz_mw_df.loc[r] = bool_parser.evaluate(enzyme_mw_dict)
            
        enz_mw_df.reset_index(inplace=True)
    
        # convert counts/fL to M (intercellular)
        # first multiplying by 10^15 (counts/fL to counts/L)
        # then dividing by Avogadro's number (counts/L to moles/L)
        enz_molar = enz_count_mean * 1e15 / 6.02e23
    
        enz_conc_df = pd.DataFrame(index=reaction2boolparser.keys(),
                                   columns=enz_count_mean.columns)
        enz_conc_df.index.name = 'bigg.reaction'
        
        for c in enz_count_mean.columns:
            enzyme_conc_dict = enz_molar[c].to_dict()
            
            for r, bool_parser in reaction2boolparser.iteritems():
                enz_conc_df.loc[r, c] = bool_parser.evaluate(enzyme_conc_dict)
        
        enz_conc_df = pd.melt(enz_conc_df.reset_index(), id_vars='bigg.reaction',
                              var_name='condition', value_name='enzyme conc [M]')
        return enz_conc_df, enz_mw_df
    
    def get_flux_data(self):
        """
            Read the flux data that was precalculated using FBA for many conditions
        """
        flux_df = pd.DataFrame.from_csv(S.ECOLI_FLUX_FNAME)
        flux_df.index.name = 'bigg.reaction'
        flux_df.reset_index(inplace=True)
        flux_df['bigg.reaction'] = flux_df['bigg.reaction'].apply(str.lower)
        flux_df = pd.melt(flux_df, id_vars='bigg.reaction',
                          var_name='condition', value_name='flux [mmol/gCDW/h]')
        return flux_df

#%%
##############################################################################

if __name__ == '__main__':
    sns.set()
    vol = Volume()
    data = vol.data_df
    E = data['enzyme conc [M]']
    r_V = data['enzyme MW [Da]'] / data['metabolite MW [Da]']
    s = data['metabolite conc [M]']
    Km = data['Km [M]']
    
    s_over_Km = r'$\frac{s}{K_M}$'
    data[s_over_Km] = s/Km
    rV_E_over_Km = r'$\frac{V_E E}{V_s K_M}$'
    data[rV_E_over_Km] = r_V * E / Km
    
    # for each metabolite:condition, take only the enzyme which
    # corresponds to the reaction with the highest flux (i.e.
    # the 'volume-limiting' reaction)
    data.sort_values('flux [mmol/gCDW/h]', inplace=True, ascending=False)
    maxflux_data = data.groupby(('bigg.metabolite', 'condition')).first().reset_index()
    
    #%%
    for hue in ['bigg.metabolite', 'bigg.reaction', 'condition']:
        pal = sns.husl_palette(len(maxflux_data[hue].unique()))
        g = sns.FacetGrid(maxflux_data, col=None, hue=hue, palette=pal,
                          ylim=(1e-3, 1e3))
        g = g.map(plt.scatter, s_over_Km, rV_E_over_Km).add_legend()
        g.ax.set_xscale('log')
        g.ax.set_yscale('log')
        g.fig.set_size_inches(12, 7)
        
        # plot the line y = x*(1+x)
        x = pd.np.logspace(-3, 3, 100)
        y = x*(1+x)
        g.ax.plot(x, y, '-')
        g.ax.set_xlim(1e-3, 1e4)
        g.ax.set_ylim(1e-3, 1e4)
        g.fig.savefig(os.path.join(S.RESULT_DIR, 'volume_%s.svg' % hue))
        
    #%%
    data['predicted_enzyme'] = 1/r_V * s * (1 + s/Km)
    data.sort_values('flux [mmol/gCDW/h]', inplace=True, ascending=False)
    maxflux_data = data.groupby(('bigg.metabolite', 'condition')).first().reset_index()
    for hue in ['bigg.metabolite', 'bigg.reaction', 'condition']:
        pal = sns.husl_palette(len(maxflux_data[hue].unique()))
        g = sns.FacetGrid(data, col=None, hue=hue, palette=pal,
                          ylim=(1e-3, 1e3))
        g = g.map(plt.scatter, 'enzyme conc [M]', 'predicted_enzyme').add_legend()
        g.ax.set_xscale('log')
        g.ax.set_yscale('log')
        g.fig.set_size_inches(12, 7)
        
        # plot the line y = x*(1+x)
        #x = pd.np.logspace(-3, 3, 100)
        #y = x*(1+x)
        #g.ax.plot(x, y, '-')
        g.ax.set_xlim(1e-9, 1e1)
        g.ax.set_ylim(1e-9, 1e1)
        g.fig.savefig(os.path.join(S.RESULT_DIR, 'volume_%s.svg' % hue))
       