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

f_chrom_size = '/users/gev/kot/stage/manta/CB4856_manta/c_elegans.WS220.fa_chrom_size.txt'
d_c_size = dico_chrom(f_chrom_size)

def manta_sv_sumsum_qual(manta_vcf, t_qual_min, t_qual_max, space, hom, d_c_size):#, t_len_sv):
     
    #extract data
    callset = allel.read_vcf(manta_vcf, fields='*')
    
    dico_chrom = {}
    s = 100000
    for chrom in np.unique(callset['variants/CHROM']):
        dico_chrom[chrom]={}
        l_qual = list(range(t_qual_min, t_qual_max + 1, space))
        l_qual.append(999)
        for t_qual in l_qual :
            if chrom=='CHROMOSOME_MtDNA':
                dico_chrom[chrom][t_qual] = np.zeros(d_c_size[chrom])
            else:
                dico_chrom[chrom][t_qual] = np.zeros(int(np.ceil(d_c_size[chrom])/s)+1)
    
    for t_qual in l_qual :
        
        for sv in np.where(callset['variants/QUAL'] >= t_qual)[0]:
            
            chromo = callset['variants/CHROM'][sv]
            if callset['variants/END'][sv]!=-1:
                
                start = callset['variants/POS'][sv]
                end = callset['variants/END'][sv]
                #print('start: ',start)
                #print('end :', end)
                if chromo=='CHROMOSOME_MtDNA':
                    for nt in range(start,end):
                        dico_chrom[chromo][t_qual][nt]+=1
                else:
                    dico_chrom[chromo][t_qual][int(np.floor(sv)/s)]+=1
                    r = list(range(start,end,s))
                    r.append(end)
                    l = []
                    for nt in range(len(r)):#start,end+1,1000):
                        ide = int(np.floor(r[nt])/s)
                        if ide not in l:
                            #print(len(dico_chrom[chromo][t_qual]))
                            dico_chrom[chromo][t_qual][ide]+=1
                            l.append(ide)
            else:
                if chromo=='CHROMOSOME_MtDNA':
                    dico_chrom[chromo][t_qual][sv]+=1
                else:
                    dico_chrom[chromo][t_qual][int(np.floor(sv)/s)]+=1
    return dico_chrom
    
    #for chrom in dico_chrom.keys():
     #   coord = np.divide()
     
    
startTime = time.time()
manta_vcf = '/users/gev/kot/stage/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv.vcf'
t_qual_min = 900
t_qual_max = 999
space = 99
f_chrom_size = '/users/gev/kot/stage/manta/CB4856_manta/c_elegans.WS220.fa_chrom_size.txt'
d_c_size = dico_chrom(f_chrom_size)   
hom = True

   
dico = manta_sv_sumsum_qual(manta_vcf, t_qual_min, t_qual_max, space, hom, d_c_size)
    
    
executionTime = (time.time() - startTime)
print('Execution time in secons : ' + str(executionTime))


df_ce = pd.read_excel('/users/gev/kot/stage/manta/CB4856_manta/cemap.xlsx')
chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
def snp_cumsum_qual(vcf_cendr, d_c_size):#, t_len_sv):
    
    df = pd.read_csv(vcf_cendr,sep='\t',names=['CHROM','POS','QUAL'])
    s = 100000
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
            print(df.iloc[snp])
            dico_chrom[ch][start]+=1
        else:
            dico_chrom[ch][int(np.floor(start)/s)]+=1

    return dico_chrom


f_chrom_size1 = '/users/gev/kot/stage/manta/cendr/c_elegans.PRJNA13758.WS276.genomic.fa_chrom_size.txt'
d_c_size1 = dico_chrom(f_chrom_size1)
vcf_cendr = '/mnt/data3/kot/cendr/WI.20200815.soft-filter_snp_pos.txt'
dic = snp_cumsum_qual(vcf_cendr, d_c_size)
    
executionTime = (time.time() - startTime)
print('Execution time in secons : ' + str(executionTime))


df_ce = pd.read_excel('/users/gev/kot/stage/manta/CB4856_manta/cemap.xlsx')
chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']


for i in dico.keys():
    if i !='CHROMOSOME_MtDNA':
        y = list(range(1,d_c_size[i]+1,1))
    else:
        y = list(range(1,d_c_size[i]+1))
    for qual in dico[i].keys():
        if qual==900:
            
        #print(np.cumsum(dico[i][qual])[-1])
            plt.plot(y,np.cumsum(dico[i][qual])/np.cumsum(dico[i][qual])[-1])
    if i !='CHROMOSOME_MtDNA':
        index = chrom.index(i) + 1
        x_ce = df_ce.loc[df_ce['chrom'] == index, 'cpos']
        y_ce = df_ce.loc[df_ce['chrom'] == index, 'genetic']/df_ce.loc[df_ce['chrom'] == index, 'genetic'].iloc[-1]
        plt.plot(x_ce,y_ce)
    plt.plot(y,np.cumsum(dic[i])/np.cumsum(dic[i])[-1])
    leg = list(dico[i].keys())
    leg=(['sv'])
    leg.append('Recombination')
    leg.append('snp')
    plt.legend(leg)
    plt.xlabel('Position pb')
    plt.ylabel('Cumulative rate %')
    title = 'Cumulative rate of SNP,recombinaison and SV\nalong the ' + i
    plt.title(title)
    plt.savefig('/users/gev/kot/stage/manta/cendr/Cumulative_rate_of_SNP_and_recombinaison_and_sv_along_the_ty' + i, dpi=400,bbox_inches="tight")
    plt.show()
    
    
    
    '''
    leg.append('ce')
    plt.legend(leg)
    plt.title(i)
    plt.show()
    '''
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