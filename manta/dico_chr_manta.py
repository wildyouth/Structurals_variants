#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 13:57:13 2022

@author: kot
"""
'''
error_file='/users/gev/kot/stage/manta/run_all_res/errror.txt'
 
with open(error_file,'r') as error_f:
    for i in error_f:
        if i[0:3]=='The':
            print(i.split(sep='/')[-4])
        

file='/users/gev/kot/stage/manta/run_all_res/all_bam_real.vcf'
 
with open(file,'r') as f:
    for i in f:
        if i[0:7]=='Warning':
            print(i)
'''      
            

           
import sys
import allel
import numpy as np

def dico_chr_manta(l_vcf_list_file, bam_path, path_manta_res, path_output):
    with open(l_vcf_list_file,'r') as l_vcf:
        c = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
        d = {}
        for ch in c:
            d[ch]=[]
        
        for vcf in l_vcf:
            callset = allel.read_vcf(vcf[:-4])
            name = vcf.split(sep='/')[-4]
            
            if callset != None:
                for chromosome in np.unique(callset['variants/CHROM']):
                    d[chromosome].append(name)
              #print(vcf[:-4])
    for chrom in d.keys():
        
        #path_manta_res = '/users/gev/kot/stage/manta/manta_founders/'
        
        path_end = '/results/variants/diploidSV_inv.vcf.gz'
        
        #path_output = '/users/gev/kot/stage/manta/manta_founders/'
        
        name_sv = path_output + chrom +'_sv_inv.txt'
        
        #bam_path = '/mnt/data3/founders_bam_real/'
        
        bam_name = path_output + chrom +'_bam_path.txt'
        
        with open(name_sv,'w') as d_chr:

            for row in d[chrom]:
                d_chr.write(path_manta_res + row + path_end + '\n')
        
        with open(bam_name,'w') as p_name:

            for row in d[chrom]:
                p_name.write(bam_path + row + '.bam\n')
        
    return

def main():
    
    #check if file is given in input
    if(len(sys.argv)<5):
        if len(sys.argv)>1:
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                print("\nYou have to give the vcf_list")
                print("USAGE: ", sys.argv[0], '<vcf_list>\n')
                print('<vcf_list>       : file that listed all vcf.gz created by run vcf_gz_list_to_txt.sh')
                print('<bam_path>       : path to the directory where all bam files are stored, all bam files should be stored in the same directory')
                print('<path_manta_res> : path to the directory where all results performed by manta are stored')
                print('<path_output>    : path to output wanted directory to store results')
            else:
                print("You have to give the vcf_list, type --help or --h to see the explanation of each file needed")
                print("USAGE: ", sys.argv[0], '<vcf_list>', '<bam_path>', '<path_manta_res>', '<path_output>')
        else:
            
            print("You have to give the vcf_list, type --help or --h to see the explanation of each file needed")
            print("USAGE: ", sys.argv[0], '<vcf_list>', '<bam_path>', '<path_manta_res>', '<path_output>')
    else:
        l_vcf_list_file = sys.argv[1]
        
        bam_path = sys.argv[2]
        
        path_manta_res = sys.argv[3]
        
        path_output = sys.argv[4]
        
        dico_chr_manta(l_vcf_list_file,bam_path,path_manta_res,path_output)
            
    return

if __name__ == '__main__':
    main()
