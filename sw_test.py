#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:42:35 2022

@author: kot
"""


import subprocess
import sys
import pathlib
# If your shell script has shebang, 
# you can omit shell=True argument.

#g1 = 'WS220'
#g2 = 'CB4856_2019'
#path_g1 = '/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa'
#path_g2 = '/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019_fasta/GCA_004526295.1_ASM452629v1_genomic.fna'

#path_output = '/users/gev/kot/stage/mauve/seed_weed_ws220_cb4856/'

def s_w(start,end, path_g1, path_g2, path_output):
    #l_r = np.linspace(start,end)
    for sw in range(start,end):
        
        sw_value = str(sw)
        g1 = pathlib.Path(path_g1).stem
        g2 = pathlib.Path(path_g2).stem
        name_output = g1 + '_' + g2 + '_sw_' + sw_value
        output_xmfa = '--output=' + path_output + name_output + '.xmfa' 
        output_tree = '--output-guide-tree=' + path_output + name_output + '.tree' 
        output_backbone = '--backbone-output=' + path_output + name_output + '.backbone' 
        
        argument = '--seed-weight=' + sw_value + ' ' + output_xmfa + ' ' + output_tree + ' ' + output_backbone  + ' ' + path_g1 + ' ' + path_g2
    
        output = subprocess.run(['/export/home1/users/gev/kot/stage/mauve/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve ' + argument], shell=True)
    
    return
def main():
    
    #check if file is given in input
    if(len(sys.argv)!=6):
        if len(sys.argv)==2:
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                print("You have to give the start and end of seed_weight testing values")
                print("USAGE: ", sys.argv[0], '<start_sw>', '<end_sw>', '<ref_gen>' , '<target_gen>', '<path_output>')
        else:
            print("You have to give the start and end of seed_weight testing values")
            print("USAGE: ", sys.argv[0], '<start_sw>', '<end_sw>', '<ref_gen>' , '<target_gen>', '<path_output>')
    else:
        
        #file_path
        start = int(sys.argv[1])
        end = int(sys.argv[2])
        path_g1=sys.argv[3]
        path_g2=sys.argv[4]
        path_output=sys.argv[5]
        s_w(start,end, path_g1, path_g2, path_output)
            
    return

if __name__ == '__main__':
    main()
