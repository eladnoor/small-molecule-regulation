# Compare Ki and KM directly

import os, sys, numpy as np, pandas as pd, pdb, matplotlib.pyplot as plt, seaborn as sns, scipy.stats as st
sns.set_style('ticks')

from statsmodels.sandbox.stats.multicomp import multipletests as padjust

plt.close('all')
plt.ion()

km = pd.read_csv('../cache/ecoli_km_bigg.csv',header = 0,index_col = 0)
ki = pd.read_csv('../cache/ecoli_ki_bigg.csv',header = 0,index_col = 0)

# Drop negative values
km = km[km['Value'] > 0]
ki = ki[ki['Value'] > 0]

# Drop duplicate values
km = km.drop_duplicates(subset = ['Value', 'EC_number', 'bigg.metabolite'])
ki = ki.drop_duplicates(subset = ['Value', 'EC_number', 'bigg.metabolite'])

# Drop mutants
ki = ki[(pd.isnull(ki['Commentary'])) |
              ((ki['Commentary'].str.find('mutant') == -1) &
               (ki['Commentary'].str.find('mutation') == -1))]
               
km = km[(pd.isnull(km['Commentary'])) |
              ((km['Commentary'].str.find('mutant') == -1) &
               (km['Commentary'].str.find('mutation') == -1))]

res = pd.DataFrame()
res['KI_Values'] = ki.groupby('bigg.metabolite')['Value'].mean()
res['KI_Number'] = ki.groupby('bigg.metabolite')['EC_number'].nunique()

res['KM_Values'] = km.groupby('bigg.metabolite')['Value'].mean()
res['KM_Number'] = km.groupby('bigg.metabolite')['EC_number'].nunique()

#res['Metabolite'] = km.groupby('bigg.metabolite')['name'].first()

# Drop rows where we don't have data for both
res = res.dropna()

# Keep only metabolites with at least 2 measurements of each
res = res[res['KI_Number'] > 1]
res = res[res['KM_Number'] > 1]


res['PValue'] = np.nan

# for each metabolite, if there is sufficient data, test
for ii in res.index:
	kid = ki[ki['bigg.metabolite'] == ii]['Value']
	kmd = km[km['bigg.metabolite'] == ii]['Value']
	
	s,p = st.mannwhitneyu( kid,kmd )
	res.at[ii,'PValue'] = p

res['QValue'] = padjust(res['PValue'],method = 'fdr_bh')[1]
res = res.sort_values('PValue')

maxval = 2*np.max( [res['KI_Values'].max(),res['KM_Values'].max()] )
minval = 0.5*np.min( [res['KI_Values'].min(),res['KM_Values'].min()] )

f,ax = plt.subplots(figsize = (8,8))
plt.loglog(res['KI_Values'],res['KM_Values'],'o')
plt.axis([minval,maxval,minval,maxval])

# Change index for readability
#res['KEGGID'] = res.index
#res.index = res['Metabolite']

# Calculate log ratio of values
res['Ratio'] = np.log2( res['KI_Values'] / res['KM_Values'] )

for ii in res.index:
	if res.at[ii,'QValue'] < 0.1:
		plt.text(res.at[ii,'KI_Values'],res.at[ii,'KM_Values'],ii)
plt.xlabel('Mean KI')
plt.ylabel('Mean KM')

diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

# The first plot is misleading, make a volcano plot
res['-Log10Q'] = -np.log10(res['QValue'])
plt.figure('Volcano')
plt.plot(res['Ratio'],res['-Log10Q'],'o')
