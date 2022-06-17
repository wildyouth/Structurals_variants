#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:31:34 2022

@author: kot
"""

'''
csv to sv to compare
'''
import pandas as pd
import numpy as np

import sys


def dico_chrom(f_chrom_size):
    l_pos = []
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
    with open(f_chrom_size,'r') as f:
        
        for line in f :
            l = line.split(sep='\t')
            l_pos.append(int(l[1]))
    pos_chrom = np.cumsum(l_pos)
    
    dico = {}
    for i in range(len(l_pos)):
        dico[chrom[i]] = pos_chrom[i]
    
    return dico

def mauve_csv_sv(file,dico):

    
    header = ['CHROM','POS','SVTYPE','END','leng']
    tab = pd.DataFrame(columns=header)
    
    df = pd.read_csv(file,usecols=['gen_index','start_lcb','end_lcb','lcb_length','sv'])
    df.drop(df[df.gen_index == 2].index, inplace=True)
    df = df.iloc[:, 1:]
    
    for i in range(len(df)):
        
        start = df['start_lcb'][i]
        end = df['end_lcb'][i]
        chrom_start_index = np.argwhere(np.array(list(dico.values()))>start)[0][0]
        
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
    
    return tab

def main():
    
    #check if file is given in input
    if(len(sys.argv)!=3):
        if len(sys.argv)>1:
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                print("You have to give the reference genome fasta file path and the csv file produced by json_to_csv.py")
                print("USAGE: ", sys.argv[0], '<path_ref_gen>', 'csv_file')
            else:
                print("You have to give the reference genome fasta file path and the csv file produced by json_to_csv.py")
                print("USAGE: ", sys.argv[0], '<path_ref_gen>', 'csv_file')
        else:
            print("You have to give the reference genome fasta file path and the csv file produced by json_to_csv.py")
            print("USAGE: ", sys.argv[0], '<path_ref_gen>', 'csv_file')
    else:
        
        path_ref_gen = sys.argv[1]
        file = sys.argv[2]
        f_chrom_size = path_ref_gen + '_chrom_size.txt'
        
        dico = dico_chrom(f_chrom_size)
        
        df_data = mauve_csv_sv(file,dico)
        
        file_name_sv = file.split('.csv')[0] + '_compare_sv.csv'
        df_data.to_csv(file_name_sv)
        
    return

if __name__ == '__main__':
    main()



