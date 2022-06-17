#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:26:00 2022

@author: kot
"""

'''
csv to sv to compare
'''
import pandas as pd
import numpy as np

def mauve_csv_sv(file):
    
    
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']


    pos_chrom = np.cumsum([15072423,15279345,13783700,17493793,20924149,17718866,13794])
    #cut -f1,2 c_elegans.WS220.dna.fa.fai > sizes.genome
    
    dico = {}
    for i in range(len(chrom)):
        dico[chrom[i]] = pos_chrom[i]
    
    header = ['CHROM','POS','SVTYPE','END','leng']
    tab = pd.DataFrame(columns=header)
    
    df = pd.read_csv(file,usecols=['gen_index','start_lcb','end_lcb','lcb_length','sv'])
    df.drop(df[df.gen_index == 2].index, inplace=True)
    df = df.iloc[:, 1:]
    
    for i in range(len(df)):
        
        start = df['start_lcb'][i]
        end = df['end_lcb'][i]
        
        chrom_start_index = np.argwhere(np.array(list(dico.values()))>start)[0][0]
        print(end)
        chrom_end_index = np.argwhere(np.array(list(dico.values()))>=end)[0][0]
        
        if chrom_start_index==chrom_end_index:
            nb_chrom = list(dico.keys())[chrom_start_index]
        else:
            nb_chrom = ''
            for j in range(chrom_start_index,chrom_end_index+1):
                nb_chrom+=list(dico.keys())[j]+'_'
            nb_chrom = nb_chrom[:-1]
        
        if chrom_start_index !=0 :
            start -= list(dico.values())[chrom_start_index-1]
        if chrom_end_index !=0 :
            end -=  list(dico.values())[chrom_end_index-1]
            
        
        tab.loc[len(tab.index)] = [nb_chrom,start,df['sv'][i],end,df['lcb_length'][i]]
    
    
    print(tab)
    
    return tab