# Compare KM and KI 

import os, numpy as np, scipy.stats as st, matplotlib.pyplot as plt
import seaborn as sns, pandas as pd, settings
import statsmodels.stats.multitest as smm

sns.set_style('ticks')

plt.ion()
plt.close('all')

# Read the KM and KI data
km = pd.read_csv(os.path.join(settings.DATA_DIR, 'ecoli_km_kegg.csv'),
                 header=0, index_col=0)
km['Type'] = 'KM'
ki = pd.read_csv('../cache/ecoli_ki_kegg.csv',header = 0, index_col = 0)
ki['Type'] = 'KI'

# Disambiguate indices, hate numerical index name
km.index = [km.ix[ii,'KEGG_ID'] + ':' + str(ii) for ii in km.index]
ki.index = [ki.ix[ii,'KEGG_ID'] + ':' + str(ii) for ii in ki.index]

# Remove any negative values
km = km[km['Value'] > 0]
ki = ki[ki['Value'] > 0]

mgroups = km.groupby('KEGG_ID')
igroups = ki.groupby('KEGG_ID')
ixmets = mgroups.groups.viewkeys() & igroups.groups.viewkeys()

# For each metabolite in ixmets, test whether there is a significant difference 
# between KM and KI. Restrict to cases where we have at least 3 samples of KM and KI.
res = pd.DataFrame( columns = ['Log2FC','MannWhitneyP','PAdj','NumKM','NumKI'] )
for met in ixmets:

    ixm = mgroups.groups[met]
    ixi = igroups.groups[met]
    
    numm = len( ixm )
    numi = len( ixi )
    
    if numm > 2 and numi > 2:
        
        u,p = st.mannwhitneyu( km.ix[ ixm,'Value' ], ki.ix[ ixi,'Value' ] )
        fc = np.log2( km.ix[ ixm,'Value' ].mean()/ki.ix[ ixi,'Value' ].mean() )
        res.ix[ met, ['Log2FC','MannWhitneyP','NumKM','NumKI'] ] = [fc,p,numm,numi]

res['PAdj'] = smm.multipletests(res['MannWhitneyP'], method='fdr_bh')[1]
res.sort_values(by = 'PAdj', ascending = True,inplace = True)

# Make a plot comparing the median KI and KM values
km_med = mgroups.median()
ki_med = igroups.median()
medplot = pd.DataFrame( {'KM':km_med.ix[ixmets,'Value'],'KI':ki_med.ix[ixmets,'Value']} )

plt.figure()
plt.loglog( medplot['KM'],medplot['KI'],'o' )
plt.xlabel('$K_M$')
plt.ylabel('$K_I$')
settings.plotdiag()

