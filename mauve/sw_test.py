#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:42:35 2022

@author: kot
"""


import subprocess
import sys
import pathlib

"""
Function to run progressiveMauve on 2 genomes reference and target
with a range of seed_weight (from start to end)

--seed-weight=<number> Use the specified seed weight for calculating initial anchors

arguement : 
    
    start : starting value of seed_weight to test
    end : ending value of seed_weight to test
    
    path_ref_gen : path of reference genome
    
    path_target_gen : path of target genome
    
    path_output : path to save progressiveMauve results
    
    path_p_mauve : path of progressiveMauve
"""

def seed_weight_run(start,end, path_ref_gen, path_target_gen, path_output, path_p_mauve):
    
    #Convert start/end to integer for loop for
    start = int(start)
    end  = int(end)
    
    for sw in range(start, end + 1): #sw = seed_weight
        
        #concert sw to string for name_output and seed_weight_arg
        sw_value = str(sw)
        
        #Get name of ref and target genome
        ref_gen = pathlib.Path(path_ref_gen).stem
        target_gen = pathlib.Path(path_target_gen).stem
        
        #Set up name of file output with ref and target name and seed_weight value tested
        name_output = ref_gen + '_' + target_gen + '_sw_' + sw_value
        
        #seed_weight value tested
        seed_weight_arg = '--seed-weight=' + sw_value
        
        #Output file : wfma, tree and backbone
        output_xmfa = '--output=' + path_output + name_output + '.xmfa' 
        output_tree = '--output-guide-tree=' + path_output + name_output + '.tree' 
        output_backbone = '--backbone-output=' + path_output + name_output + '.backbone' 
        
        #list of all argument to perform progressiveMauve
        #progressiveMauve --seed_weight=sw --output=xmfa_file --output-guide-tree=guide_tree_file --backbone-output=backbone_file ref_gen target_gen
        
        l_arg = [path_p_mauve, seed_weight_arg ,output_xmfa, output_tree, output_backbone, path_ref_gen, path_target_gen]
        
        #join with space to perform in subprocess.run
        arg = ' '.join(l_arg)
        
        #run progressiveMauve with given parameters
        subprocess.run([arg], shell=True)
        
    return

def main():
    
    #if all 6 arguements are not given
    if(len(sys.argv)!=7):
        
        if len(sys.argv)==2:
            
            #print ^assistant to help to perform this script
            if sys.argv[1]=='--help' or sys.argv[1]=='--h':
                
                print("You have to give the start and end of seed_weight testing values\n")
                print("USAGE: ", sys.argv[0], '<start_sw>', '<end_sw>', '<ref_gen>' , '<target_gen>', '<path_output>', '<path_progressive_mauve>\n')
                print("start_sw : starting value of seed_weight to test")
                print("end_sw : ending value of seed_weight to test")
                print("ref_gen : path of reference genome")
                print("target_gen : path of target genome")
                print("path_output : path to save progressiveMauve results")
                print("path_progressive_mauve : path of progressiveMauve")
                
        else:
            
            print("You have to give the start and end of seed_weight testing values")
            print("USAGE: ", sys.argv[0], '<start_sw>', '<end_sw>', '<ref_gen>' , '<target_gen>', '<path_output>', 'path_progressive_mauve')
            print("USAGE: ", sys.argv[0], '--h to have arguments details')
    
    #all 6 arguemnts are given
    else:
        
        start, end, path_ref_gen, path_target_gen, path_output, path_p_mauve = sys.argv[1:]
        
        #run progressiveMauve
        seed_weight_run(start,end, path_ref_gen, path_target_gen, path_output, path_p_mauve)
            
    return

if __name__ == '__main__':
    main()