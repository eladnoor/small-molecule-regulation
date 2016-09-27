# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 19:49:04 2016

@author: noore
"""
import pandas as pd
import settings as S
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from bool_parser import BoolParser

def get_metabolite_data():
    """
        get the metabolomics data from Gerosa et al. 2015
    """
    _df = pd.DataFrame.from_csv(S.ECOLI_METAB_FNAME)
    _df.index.name = 'bigg.metabolite'
    met_conc_mean = _df.iloc[:, 1:9] * 1e-3 # convert mM to M
    met_conc_std = _df.iloc[:, 10:]  * 1e-3 # convert mM to M
    
    met_conc_mean.columns = [c[:-7].lower() for c in met_conc_mean.columns]
    met_conc_std.columns = [c[:-6].lower() for c in met_conc_std.columns]
    
    return met_conc_mean

def get_km_data():
    km_df = S.read_cache('km')
    km_df = km_df[km_df['Organism'] == 'Escherichia coli']
    km_df = km_df[km_df['KM_Value'] != -999]
    km_df = km_df[['EC_number', 'KM_Value', 'bigg.metabolite']]
    
    # join with ec2bigg table to map values to bigg.reactions
    ec2bigg = pd.DataFrame.from_csv(S.BIGG2EC_FNAME)
    km_df = pd.merge(km_df, ec2bigg, on='EC_number')
    return km_df

def get_enzyme_data():
    """
        get the enzyme data from Schmidt et al. 2015
    """
    _df = pd.DataFrame.from_csv(S.ECOLI_PROT_FNAME)
    _df = _df[~pd.isnull(_df['Bnumber'])]
    _df.set_index('Bnumber', inplace=True)
    enz_mw         = _df['Molecular weight (Da)']
    enz_count_mean = _df.iloc[:, 7:29] # counts per 10^-15 liter
    enz_count_cv   = _df.iloc[:, 30:52]
    enz_count_mean.columns = [c[:-7].lower() for c in enz_count_mean.columns]
    enz_count_cv.columns = [c[:-5].lower() for c in enz_count_cv.columns]
    
    enz_count_mean = enz_count_mean.apply(pd.to_numeric).fillna(0.0)
    
    # for each reaction, get a single value for the enzyme concentration
    # in [M] and the molecular weight of the complex
    
    model = create_cobra_model_from_sbml_file(S.ECOLI_SBML_FNAME)
    convert_to_irreversible(model)
    
    reaction2boolparser = {}
    for r in model.reactions:
        if r.gene_reaction_rule == '':
            continue
        if r.gene_reaction_rule.find('s0001') != -1:
            continue
        reaction2boolparser[r.id.lower()] = BoolParser(r.gene_reaction_rule)
    
    enz_conc_df = pd.DataFrame(index=reaction2boolparser.keys(),
                               columns=enz_count_mean.columns)
    enz_mass_df = pd.DataFrame(index=reaction2boolparser.keys(),
                               columns=enz_count_mean.columns)
    enz_conc_df.index.name = 'bigg.reaction'
    enz_mass_df.index.name = 'bigg.reaction'

    
    rho = 1100 # average cell density gr/L
    DW_fraction = 0.3 # fraction of DW of cells
    enz_molar = enz_count_mean * 1e15 / 6.02e23 / rho # convert counts/fL to M (intercellular)
    
    for c in enz_count_mean.columns:
        enz_gr_per_gCDW = enz_molar[c] * enz_mw / DW_fraction # convert M to gr/gCDW
        enzyme_conc_dict = enz_molar[c].to_dict()
        enzyme_mass_dict = enz_gr_per_gCDW.to_dict()
        
        for r, bool_parser in reaction2boolparser.iteritems():
            enz_conc_df.loc[r, c] = bool_parser.evaluate(enzyme_conc_dict)
            enz_mass_df.loc[r, c] = bool_parser.evaluate(enzyme_mass_dict)
    
    return enz_conc_df, enz_mass_df


##############################################################################
met_conc_df = get_metabolite_data()
km_df = get_km_data()
enz_conc_df, enz_mass_df = get_enzyme_data()

# get the flux from FBA
flux_df = pd.DataFrame.from_csv(S.ECOLI_FLUX_FNAME)
flux_df.index.name = 'bigg.reaction'
flux_df.reset_index(inplace=True)
flux_df['bigg.reaction'] = flux_df['bigg.reaction'].apply(str.lower)

#%%
enz_conc_df = pd.melt(enz_conc_df.reset_index(), id_vars='bigg.reaction',
                      var_name='condition', value_name='enzyme conc [M]')
enz_mass_df = pd.melt(enz_mass_df.reset_index(), id_vars='bigg.reaction',
                      var_name='condition', value_name='enzyme mass [gr/gCDW]')
met_conc_df = pd.melt(met_conc_df.reset_index(), id_vars='bigg.metabolite',
                      var_name='condition', value_name='metabolite conc [M]')
flux_df = pd.melt(flux_df, id_vars='bigg.reaction',
                      var_name='condition', value_name='flux [mmol/gCDW/h]')

#%%
data_df = pd.merge(km_df, met_conc_df, on='bigg.metabolite')
data_df = pd.merge(data_df, enz_conc_df, on=['bigg.reaction', 'condition'])
data_df = pd.merge(data_df, enz_mass_df, on=['bigg.reaction', 'condition'])
data_df = pd.merge(data_df, flux_df, on=['bigg.reaction', 'condition'])

# keep only rows with non-zero flux
data_df = data_df[data_df['flux [mmol/gCDW/h]'] != 0]
# keep only rows with non-zero enzyme concentration
data_df = data_df[data_df['enzyme conc [M]'] > 0]


# keep only the intersection of all conditions
# for which we have metabolites, proteins and flux
#cond = met_conc_df.columns & enz_conc_df.columns & flux_df.columns
