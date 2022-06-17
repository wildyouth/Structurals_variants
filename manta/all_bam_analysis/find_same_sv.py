#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:19:02 2022

@author: kot
"""

import pandas as pd

def is_unique(s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()


def df_to_keep (df_all, chrom,start,sv_type,end,sv_len,start_space,end_space,sv_len_range):
    
    if sv_type != 'BND':
        
        df = df_all.loc[(df_all['CHROM'] == chrom) & 
                            (df_all['POS'] < start + start_space) & (df_all['POS'] > start - start_space) & 
                            (df_all['END'] < end + end_space) & (df_all['END'] > end - end_space) & 
                            (df_all['SVTYPE'] == sv_type) & 
                            (df_all['SVLEN'] > sv_len* 1-sv_len_range) &
                            (df_all['SVLEN'] < sv_len* 1+sv_len_range) ]
    
    else:
        df = df_all.loc[(df_all['CHROM'] == chrom) & 
                            (df_all['POS'] < start + start_space) & (df_all['POS'] > start - start_space) & 
                            (df_all['END'] < end + end_space) & (df_all['END'] > end - end_space) & 
                            (df_all['SVTYPE'] == sv_type) & 
                            (df_all['SVLEN'] == -1) ]
        
    return df

file = '/users/gev/kot/stage/results_gt/CHROMOSOME_all_vcf_filter_true.csv'

df_all = pd.read_csv(file)

manta_cb_4856_file = '/users/gev/kot/stage/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv_compare_sv_qual_900_hom_False.csv'

df_cb4856 = pd.read_csv(manta_cb_4856_file)

#col = ,CHROM,POS,ID,END,FILTER_PASS,SVTYPE,SVLEN

header = ['CHROM','POS','SVTYPE','END']
col = ['CHROM','POS','ID','END','FILTER_PASS','SVTYPE','SVLEN']
df  = pd.DataFrame(columns=col)
df_check = pd.DataFrame(columns=col)
for idx, row in df_cb4856.iterrows():
    
    chrom = row['CHROM']
    start = row['POS']
    sv_type = row['SVTYPE']
    end = row['END']
    sv_len = row['SVLEN']
    if sv_type == 'DEL':
        sv_len = -sv_len
    start_space = 200
    end_space = 200
    sv_len_range = 0.2
    
    add_df = df_to_keep(df_all, chrom,start,sv_type,end,sv_len,start_space,end_space,sv_len_range)
    print(add_df)
    if (len(add_df.index) > 1) :
        a = is_unique(add_df['POS']) & is_unique(add_df['END']) & is_unique(add_df['SVLEN'])
        if not a :
            start_space = 100
            end_space = 100
            sv_len_range = 0.1
            add_df = df_to_keep(df_all, chrom,start,sv_type,end,sv_len,start_space,end_space,sv_len_range)
            if (len(add_df.index) > 1) :
                a = is_unique(add_df['POS']) & is_unique(add_df['END']) & is_unique(add_df['SVLEN'])
                if not a :
                    add_df = add_df.loc[:, add_df.columns != 'Unnamed: 0']
    
                    df_check = pd.concat([df_check,add_df])
                    print(add_df)
        else:
            add_df=add_df.head(1)
    add_df = add_df.loc[:, add_df.columns != 'Unnamed: 0']
    
    df = pd.concat([df,add_df])
    
df.to_csv('/users/gev/kot/stage/'+'CHROMOSOME_all_vcf_same_sv_all.csv')
df_check.to_csv('/users/gev/kot/stage/'+'CHROMOSOME_all_check_vcf_same_sv_all.csv')