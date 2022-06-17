#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 11:41:25 2022

@author: kot
"""

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
#pip3 install scikit-allel
import allel
import time
#print(manta_vcf_to_csv(manta_vcf, 200, False))#, -1))
'''
return df of diploidSV_inv.vcf to save in csv after
'''
'''
def dico_chrom(f_chrom_size):
    l_pos = []
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
    with open(f_chrom_size,'r') as f:
        
        for line in f :
            l = line.split(sep='\t')
            l_pos.append(int(l[1]))
    #pos_chrom = np.cumsum(l_pos)
    
    dico = {}
    for i in range(len(l_pos)):
        dico[chrom[i]] = l_pos[i]
    
    return dico

f_chrom_size = '/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa_chrom_size.txt'
d_c_size = dico_chrom(f_chrom_size)
'''


def manta_select_large_sv(manta_vcf, size_min,t_qual, hom):
    #extract data
    callset = allel.read_vcf(manta_vcf, fields='*')
    
    #columns
    header = ['CHROM','POS','SVTYPE','END']
    
    df  = pd.DataFrame(columns=header)
    
    gt = allel.GenotypeArray(callset['calldata/GT'])
    
    #threshold selct only rows with qual greater than t_qual
    thres_qual = np.where(callset['variants/QUAL'] > t_qual)
    leng = callset['variants/END']-callset['variants/POS']
    t_size = np.where(leng>size_min)
    
    #check if user wants to get ride of non homologus
    if hom == True :
        l_homo = np.where(gt.is_hom() == hom)
        param = np.intersect1d(l_homo[0], thres_qual[0],t_size)
    else:
        param = np.intersect1d(thres_qual[0],t_size)
    
    #add data to df
    for col in header :
        df[col]=np.take(callset['variants/' + col], param)
    
    #df['leng'] = df['END'] - df['POS']
    
    #df = df.assign(leng=(df['POS'] != df['END']))
    #df = df.assign(leng=(df['END'] - df['POS'] >= t_len_sv) or (df['END'] == -1))
    
    #df = df.loc[df['leng'] == True]
    
    #df = df.iloc[:, 0:4]
    print(df)
    return df

manta_vcf = '/users/gev/kot/stage/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv.vcf'
t_qual = 900
hom=False
size_min = 1000000
d = manta_select_large_sv(manta_vcf, size_min, t_qual, hom)
d.to_csv('/users/gev/kot/stage/manta/sv_cb4856_manta_large_sv_900.csv')
'''

def main():
    
    param = [0, False]
    
    if(len(sys.argv)!=1):
        #check if file is given in input
        if(len(sys.argv) == 2):
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                print("You have to give manta vcf file")
                print("USAGE: ", sys.argv[0], '[options]', '<manta vcf file>')
                print('Options:\n')
                print('\t --thres_qual = <number> Threshold value of QUAL param of vcf, default value [0]')
                print('\t --hom = Bool True/False Get rid (True) or not (False) of non homologus read, default value [False]')
            else:
                manta_vcf = sys.argv[1]
            
        elif (len(sys.argv) >4):
            
            print("ERROR too many arguments")
        
        else:
            for arg in sys.argv[:-1] :
                option = arg.split(sep='=')
                if option[0]=='--thres_qual':
                    param[0]=option[1]
                elif option[0]=='--hom':
                    param[1]=option[1]
            
            manta_vcf = sys.argv[-1]
            
        df_data = manta_vcf_to_csv(manta_vcf, float(param[0]), param[1])
        
        file_name_sv = manta_vcf.split('.vcf')[0] + '_compare_sv_'+'qual_'+str(param[0])+'_hom_'+str(param[1])+'.csv'
        df_data.to_csv(file_name_sv)
    else:
        print("You have to give manta vcf file")
        print("USAGE: ", sys.argv[0], '[options]', '<manta vcf file>')
        
    return

if __name__ == '__main__':
    main()
'''