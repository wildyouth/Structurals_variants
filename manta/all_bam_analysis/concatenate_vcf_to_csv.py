#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 17:48:35 2022

@author: kot
"""
import re
import os
from functools import reduce
import pandas as pd
import allel

path_all_chr_with_vcf = '/users/gev/kot/stage/results_gt/CHROMOSOME_MtDNA/'

def get_files(path_all_chr_with_vcf):
    
    #get all csv files in given path in a list 
    #for root, dirs, files in os.walk(path):
    #    l_files = files.copy()
    l_summary_files = []
    l_summary_files += [file for file in os.listdir(path_all_chr_with_vcf) if file.endswith('.vcf')]

    if len(l_summary_files)==0:
        raise Exception("there is 0 vcf file in given directory")
    return l_summary_files
all_df = pd.DataFrame()
for i in ['I','II','III','IV','V','X','MtDNA']:
    
    #all_df = pd.DataFrame()
    
    path_all_chr_with_vcf = '/users/gev/kot/stage/results_gt/CHROMOSOME_' + i + '/'

    l_summary_files = get_files(path_all_chr_with_vcf)
    print(l_summary_files)
    
    

    for vcf in l_summary_files :
        
        callset = allel.read_vcf(path_all_chr_with_vcf + vcf, fields='*')#['samples','CHROM','POS','ID','END','FILTER_PASS','SVTYPE','SVLEN','GT',])
        df = allel.vcf_to_dataframe(path_all_chr_with_vcf + vcf, fields='*')#['samples','CHROM','POS','ID','END','FILTER_PASS','SVTYPE','SVLEN','GT',])
        add_df = df.copy()
        all_df = all_df.append(add_df, ignore_index=True)
        
    #print(all_df)
    #final_df =  reduce(lambda  left,right: pd.merge(left,right,on=['CHROM'],how='outer'), all_df)
    #print(final_df)
all_df.to_csv(path_all_chr_with_vcf+'CHROMOSOME_all_vcf_all_data.csv')


file = '/users/gev/kot/stage/results_gt/CHROMOSOME_MtDNA/000000001-000013794.vcf'
callset = allel.read_vcf(file, fields='*')
samples = callset['samples']
for i in 



