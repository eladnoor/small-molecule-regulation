# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 14:39:33 2016

@author: noore
"""

from bioservices.kegg import KEGG
import pandas as pd
import re, os
import settings

k = KEGG()
k.organism = 'eco'

def get_ec_df():
    eco2ec = k.link('eco', 'ec')
    eco2ec = filter(lambda x: len(x) == 2, map(lambda l : l.split('\t'), eco2ec.split('\n')))
    eco2ec = pd.DataFrame(eco2ec)
    eco2ec[0] = eco2ec[0].apply(lambda s: s[3:])
    eco2ec[1] = eco2ec[1].apply(lambda s: s[4:])
    eco2ec.rename(columns={0: 'EC_number', 1: 'b_number'}, inplace=True)
    eco2ec.set_index('EC_number', inplace=True)
    eco2ec.index.name = 'EC_number'
    return eco2ec

def get_kegg_df():
    if not os.path.exists(settings.KEGG2CHEBI_FNAME):
        cid2chebi = k.list('cpd')
        cid2chebi = filter(lambda x: len(x) == 2, map(lambda l : l.split('\t'), cid2chebi.split('\n')))
        cid2chebi = pd.DataFrame(cid2chebi)
        cid2chebi[0] = cid2chebi[0].apply(lambda s: s[4:])
        cid2chebi[1] = cid2chebi[1].apply(lambda s: s.split(';')[0])
        cid2chebi.rename(columns={0: 'KEGG_ID', 1: 'Name'}, inplace=True)
        cid2chebi.set_index('KEGG_ID', inplace=True)
        cid2chebi.index.name = 'KEGG_ID'
        
        cid2chebi['ChEBI'] = None
        for cid in cid2chebi.index[0:10]:
            ChEBI = re.findall('ChEBI: ([\d\s]+)\n', k.get(cid))
            if len(ChEBI) == 0:
                print 'Cannot find a ChEBI for %s' % cid
            elif len(ChEBI) > 1:
                print 'Error parsing compound %s' % cid
            else:
                cid2chebi.at[cid, 'ChEBI'] = ChEBI[0]
        
        cid2chebi.to_csv(settings.KEGG2CHEBI_FNAME)
    
    kegg_df = pd.DataFrame.from_csv(settings.KEGG2CHEBI_FNAME,
                                    header=0, index_col=None)
    
    # remove compounds that have no ChEBI:
    kegg_df = kegg_df[~kegg_df['ChEBI'].isnull()]
    kegg_df.set_index('KEGG_ID', inplace=True)
    
    # split compounds with more than one ChEBI to several rows:
    tmp_chebi_df = pd.DataFrame(index=kegg_df.index, 
                                data=kegg_df['ChEBI'].apply(str.split).tolist())
    kegg_df = kegg_df.join(tmp_chebi_df)
    kegg_df = kegg_df.drop('ChEBI', axis=1)
    kegg_df['KEGG_ID'] = kegg_df.index
    kegg_df = pd.melt(kegg_df, id_vars=['KEGG_ID', 'name'], value_name='ChEBI')
    kegg_df = kegg_df.drop('variable', axis=1)
    kegg_df = kegg_df[~kegg_df['ChEBI'].isnull()]
    kegg_df.sort_values('KEGG_ID', inplace=True)
    return kegg_df