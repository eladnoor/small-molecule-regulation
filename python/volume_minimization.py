# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 19:49:04 2016

@author: noore
"""
import json
import pandas as pd
import settings as S
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from bool_parser import BoolParser


# load the metabolomics data from Gerosa et al. 2015
_df = pd.DataFrame.from_csv(S.ECOLI_METAB_FNAME)
_df.index.name = 'bigg.metabolite'
met_conc_mean = _df.iloc[:, 1:9]
met_conc_std = _df.iloc[:, 10:]

met_conc_mean.columns = [c[:-7].lower() for c in met_conc_mean.columns]
met_conc_std.columns = [c[:-6].lower() for c in met_conc_std.columns]


# load the proteomic data from Schmidt et al. 2015
_df = pd.DataFrame.from_csv(S.ECOLI_PROT_FNAME)
_df = _df[~pd.isnull(_df['Bnumber'])]
_df.set_index('Bnumber', inplace=True)
enz_mw         = _df['Molecular weight (Da)']
enz_count_mean = _df.iloc[:, 7:29] # counts per 10^-15 liter
enz_count_cv   = _df.iloc[:, 30:52]
enz_count_mean.columns = [c[:-7].lower() for c in enz_count_mean.columns]
enz_count_cv.columns = [c[:-5].lower() for c in enz_count_cv.columns]

enz_count_mean = enz_count_mean.apply(pd.to_numeric).fillna(0.0)

# merge with Km data
km = S.read_cache('km')
km = km[km['Organism'] == 'Escherichia coli']
km = km[km['KM_Value'] != -999]

data = km[['EC_number', 'KM_Value', 'bigg.metabolite']].join(met_conc_mean, on='bigg.metabolite', how='inner')
data.set_index('bigg.metabolite', inplace=True)

# now we need to map the enzyme counts to EC numbers and join with the
# metabolite+Km data table

# for that, we need to use the E. coli model to map between bnumbers and 
# reactions and their EC-numbers, we also need to make sure we only take
# reactions with positive flux (in the irreversible model)

flux_df = pd.DataFrame.from_csv(S.ECOLI_FLUX_FNAME)
flux_df.index.name = 'bigg.reaction'

# common conditions:
cond = met_conc_mean.columns & enz_count_mean.columns & flux_df.columns
met_conc_mean = met_conc_mean[cond]
met_conc_std = met_conc_std[cond]
enz_count_mean = enz_count_mean[cond]
enz_count_cv = enz_count_cv[cond]
flux_df = flux_df[cond]

model = create_cobra_model_from_sbml_file(S.ECOLI_SBML_FNAME)
convert_to_irreversible(model)

gene_df = pd.DataFrame(data=[(r.id, r.gene_reaction_rule) for r in model.reactions],
                       columns=('bigg.reaction', 'gene_reaction_rule'))
gene_df.set_index('bigg.reaction', inplace=True)

#%% calculate the total mass of enzymes associated with each reaction in each condition
for c in cond:
    gene_df[c] = 0
    
    rho = 1100 # average cell density gr/L
    DW_fraction = 0.3 # fraction of DW of cells
    enz_gr_per_gCDW = enz_count_mean[c] * enz_mw * (1e15 / 6.02e23 / rho / DW_fraction)
    enzyme_mass_dict = enz_gr_per_gCDW.to_dict()
    
    for r in gene_df.index:
        if gene_df.loc[r, 'gene_reaction_rule'] == '':
            continue
        if gene_df.loc[r, 'gene_reaction_rule'].find('s0001') != -1:
            continue
        
        bool_parser = BoolParser(gene_df.loc[r, 'gene_reaction_rule'])
        gene_df.loc[r, c] = bool_parser.evaluate(enzyme_mass_dict)

gene_df.melt()