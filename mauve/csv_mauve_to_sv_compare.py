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

"""
Function to extract chromosome size in dictionnary

arguement : 
    
    f_chrom_size : file with chromosomes size of reference genome
"""

def dico_chrom(f_chrom_size):
    
    
    #store chromosome size in a list
    l_pos = []
    
    #name of chromosomes in C. elegans
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
    
    #open file
    with open(f_chrom_size,'r') as f:
        
        for line in f :
            l = line.split(sep='\t')
            
            #add chromosome size
            l_pos.append(int(l[1]))
    
    #progressiveMauve concatenate all chromosome coordinate
    #To give each lcbs the chromosome associated to them we need the coordinate of reference genome concatenated
    pos_chrom = np.cumsum(l_pos)
    
    dico_concanated_chrom_coord = {}
    
    
    for chrom_concatenated_coordinate in range(len(l_pos)):
        dico_concanated_chrom_coord[chrom[chrom_concatenated_coordinate]] = pos_chrom[chrom_concatenated_coordinate]
    
    return dico_concanated_chrom_coord

def mauve_csv_sv(csv_file,dico_concanated_chrom_coord):
    
    '''
    keep only wanted data
    CHROM = Chromosome in reference genome of SV
    POS = SV starting postion
    END = SV ending position
    SVTYPE = type of SV
    leng = SV length
    '''
    
    header = ['CHROM','POS','SVTYPE','END','leng']
    new_df_sv = pd.DataFrame(columns=header)
    
    #keep only wanted columns
    df_file = pd.read_csv(csv_file,usecols=['gen_index','start_lcb','end_lcb','lcb_length','sv'])
    
    #exclude target genome
    df_file.drop(df_file[df_file.gen_index == 2].index, inplace=True)
    #exlude gen_index column
    df_file = df_file.iloc[:, 1:]
    
    
    for sv in range(len(df_file)):
        
        start = df_file['start_lcb'][sv]
        end = df_file['end_lcb'][sv]
        
        #give each sv their chromosome tag
        chrom_start_index = np.argwhere(np.array(list(dico_concanated_chrom_coord.values()))>start)[0][0]
        
        chrom_end_index = np.argwhere(np.array(list(dico_concanated_chrom_coord.values()))>=end)[0][0]
        
        #check if progressiveauve find sv overlaping 2 chromosomes
        if chrom_start_index==chrom_end_index:
            nb_chrom = list(dico_concanated_chrom_coord.keys())[chrom_start_index]
        else:
            nb_chrom = ''
            for j in range(chrom_start_index,chrom_end_index+1):
                nb_chrom+=list(dico_concanated_chrom_coord.keys())[j]+'_'
            nb_chrom = nb_chrom[:-1]
        
        #give the chromosome coordinate without being concanated
        if chrom_start_index !=0 :
            start -= list(dico_concanated_chrom_coord.values())[chrom_start_index-1]
        if chrom_end_index !=0 :
            end -=  list(dico_concanated_chrom_coord.values())[chrom_end_index-1]
            
        
        new_df_sv.loc[len(new_df_sv.index)] = [nb_chrom,start,df_file['sv'][sv],end,df_file['lcb_length'][sv]]
    
    return new_df_sv

def main():
    
    #check if file is given in input
    if(len(sys.argv)!=3):
        if len(sys.argv)>1:
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                print("You have the csv file produced by json_to_csv.py and file with chromosomes size of reference genome")
                print("USAGE: ", sys.argv[0], '<csv_file>', '<path_chrom_size_ref_gen>')
                print("csv_file : csv file from json_to_csv.py")
                print("path_chrom_size_ref_gen : path of txt file with chromosomes size of reference genome")
            else:
                print("You have the csv file produced by json_to_csv.py and file with chromosomes size of reference genome")
                print("USAGE: ", sys.argv[0], '<csv_file>', '<path_chrom_size_ref_gen>')
        else:
            print("You have the csv file produced by json_to_csv.py and file with chromosomes size of reference genome")
            print("USAGE: ", sys.argv[0], '<csv_file>', '<path_chrom_size_ref_gen>')
    else:
        
        csv_file = sys.argv[1]
        f_chrom_size = sys.argv[2]
        
        dico = dico_chrom(f_chrom_size)
        
        df_data = mauve_csv_sv(csv_file,dico)
        
        file_name_sv = csv_file.split('.csv')[0] + '_compare_sv.csv'
        df_data.to_csv(file_name_sv)
        
    return

if __name__ == '__main__':
    main()



