#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:43:56 2022

@author: kot
"""

import pandas as pd
import numpy as np 
import time
import matplotlib.pyplot as plt

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

f_chrom_size = '/users/gev/kot/stage/manta/cendr/c_elegans.PRJNA13758.WS276.genomic.fa_chrom_size.txt'
d_c_size = dico_chrom(f_chrom_size)


def snp_cumsum_qual(vcf_cendr, d_c_size):#, t_len_sv):
    
    df = pd.read_csv(vcf_cendr,sep='\t',names=['CHROM','POS','QUAL'])
    s = 2500
    dico_chrom = {}

    for chrom in df.CHROM.unique():
        ch = 'CHROMOSOME_' + chrom

        if ch=='CHROMOSOME_MtDNA':

            dico_chrom[ch] = np.zeros(d_c_size[ch])

        else:
            dico_chrom[ch] = np.zeros(int(np.ceil(d_c_size[ch])/s)+1)
    for snp in range(len(df.CHROM)):
        
        start = df.POS[snp]
        ch = 'CHROMOSOME_' + df.CHROM[snp]
    
        if ch=='CHROMOSOME_MtDNA':
            #print(df.iloc[snp])
            dico_chrom[ch][start]+=1
        else:
            dico_chrom[ch][int(np.floor(start)/s)]+=1

    return dico_chrom
    
    #for chrom in dico_chrom.keys():
     #   coord = np.divide()
     
     
'''
df = df[df['QUAL'] != '.']
np.min(df.iloc[:,2:].values)
Out[66]: 30.0

np.max(df.iloc[:,2:].values)
Out[67]: 43912300.0
'''
     
    
startTime = time.time()

f_chrom_size = '/users/gev/kot/stage/manta/cendr/c_elegans.PRJNA13758.WS276.genomic.fa_chrom_size.txt'
d_c_size = dico_chrom(f_chrom_size)
vcf_cendr = '/mnt/data3/kot/cendr/WI.20200815.soft-filter_snp_pos.txt'
dico = snp_cumsum_qual(vcf_cendr, d_c_size)

    
executionTime = (time.time() - startTime)
print('Execution time in secons : ' + str(executionTime))


df_ce = pd.read_excel('/users/gev/kot/stage/manta/CB4856_manta/cemap.xlsx')
chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']




for i in dico.keys():
    if i !='CHROMOSOME_MtDNA':
        y = list(range(1,d_c_size[i]+1,2500))
    else:
        y = list(range(1,d_c_size[i]+1))
    plt.plot(y,np.cumsum(dico[i])/np.cumsum(dico[i])[-1])
    if i !='CHROMOSOME_MtDNA':
        index = chrom.index(i) + 1
        x_ce = df_ce.loc[df_ce['chrom'] == index, 'cpos']
        y_ce = df_ce.loc[df_ce['chrom'] == index, 'genetic']/df_ce.loc[df_ce['chrom'] == index, 'genetic'].iloc[-1]
        plt.plot(x_ce,y_ce)
    leg = ['SNP']
    leg.append('Recombination')
    plt.legend(leg)
    plt.xlabel('Position pb')
    plt.ylabel('Cumulative proportion %')
    title = 'Cumulative proportion of SNP and recombination\nalong the ' + i
    plt.title(title)
    plt.savefig('/users/gev/kot/stage/manta/cendr/Cumulative_rate_of_SNP_and_recombination_along_the_' + i, dpi=400)
    plt.show()
    
'''   
    
    return df

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
    main()'''