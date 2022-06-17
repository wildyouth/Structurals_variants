#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:45:37 2022

@author: kot
"""

import numpy as np
import pandas as pd
import sys
#pip3 install scikit-allel
import allel

'''
return df of diploidSV_inv.vcf to save in csv after
'''
def manta_vcf_to_csv(manta_vcf, t_qual, hom):
    
    #extract data
    callset = allel.read_vcf(manta_vcf, fields='*')
    
    #columns
    header = ['CHROM','POS','ID','END','FILTER_PASS','SVTYPE','SVLEN']
    
    df  = pd.DataFrame(columns=header)
    
    gt = allel.GenotypeArray(callset['calldata/GT'])
    
    #threshold selct only rows with qual greater than t_qual
    thres_qual = np.where(callset['variants/QUAL'] > t_qual)
    
    #check if user wants to get ride of non homologus
    if hom == True :
        l_homo = np.where(gt.is_hom() == hom)
        param = np.intersect1d(l_homo[0], thres_qual[0])
    else:
        param = np.copy(thres_qual[0])
    
    #add data to df
    for col in header :
        df[col]=np.take(callset['variants/' + col], param)

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
                df_data = manta_vcf_to_csv(manta_vcf, float(param[0]), param[1])
            
                file_name_sv = manta_vcf.split('.vcf')[0] + '_compare_sv_'+'qual_'+str(param[0])+'_hom_'+str(param[1])+'.csv'
                df_data.to_csv(file_name_sv)
            
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