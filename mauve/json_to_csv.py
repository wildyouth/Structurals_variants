#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 11:05:10 2022

@author: kot
"""


import json
import os
import pandas as pd
import sys
import numpy as np
"""
check if the input file is empty
Used in read_json function
    Argument :
        file : json_file path

    Return  True if file is empty
            False is the file is not empty 
""" 
def file_is_empty(file):
    #return True is empty and False if not
    return os.stat(file).st_size==0
"""
Extract data from json file

    Argument :
        file : json_file path
    
    return  list1 of list2 of dictionary
            list -> list of all part
            list2 -> list of all dictonnary
"""
def read_json(file): 
    
    #check file if empty
    if file_is_empty(file):
        data = {}
    else:
        #open file and load data
        with open(file,'r') as json_file:
            data = json.load(json_file)
    return data    

def CheckForMore(list1, val): 
    return(all(x >= val for x in list1))  

"""
Extract data of list from json_file
    
    Argument : 
        list_d : read_json output, list of list of dico
    
    return tab with all data cf line 60 + another column with relative index of lcb
"""
def extract_json(list_d):
    
    #Create df
    #column names
    header = ['genome','gen_index','start_lcb','end_lcb','lcb_length','nb_lcb','strand','nb_gap','l_gap','average_length_gap']
    
    tab = pd.DataFrame(columns=header)
    
    #count nb_lcb
    nb_lcb = 0
    
    #iterate in list_d to parse each lcb
    for l_lcb in list_d:
        l_length = []
        for gen in l_lcb:
            start = gen['start']        #position of 1st nt of the lcb
            end = gen['end']            #position of last nt of the lcb
            length = end-start
            l_length.append(length)
        
        #only keep LCBs longer than 50bp 
        if CheckForMore(l_length, 50):
            #if len(l_lcb)==2 it means that the lcb is in both genome
            if len(l_lcb)==2:
                
                #add +1 to counter
                nb_lcb += 1
                
                l_gen = [[],[]]
                #extract data for each genome
                n_gen = 0
                
                l_length = []
                for gen in l_lcb:
                    
                    genome = gen['name']        #Name of genome
                    gen_index = gen['lcb_idx']  #genome_index 1=reference
                    start = gen['start']        #position of 1st nt of the lcb
                    end = gen['end']            #position of last nt of the lcb
                    length = end-start          #length of lcb
                    strand  = gen['strand']     #inversion or not
                    nb_gap = len(gen['gaps'])   #nb of gaps 
                    l_length.append(length)
                    #list to store all gap's length
                    l_gap = []
                    if nb_gap !=0:
                        
                        #iterate for each gap
                        for gap in (gen['gaps']):
                            
                            #add in l_gap gap's length
                            l_gap.append(gap['end']-gap['start'])
                        
                        #average gap_length
                        a_length_gap = sum(l_gap)/nb_gap
                    
                    else:
                        #average gap_length
                        a_length_gap = 0
                    
                    l_gen[n_gen] = [genome,gen_index,start,end,length,nb_lcb,strand,nb_gap,l_gap,a_length_gap]
                    n_gen+=1
                
                if 0 in l_length :
                    
                    start = np.nan              #position of 1st nt of the lcb
                    end = np.nan                #position of last nt of the lcb
                    length = 0                  #length of lcb
                    nb_lcb = 1
                    strand  = np.nan            #inversion or not
                    nb_gap = 0                #nb of gaps 
                    
                    
                    for gen_data in l_gen :
                    #add all data in df  
                        gen_data[2:8]=[start,end,length,nb_lcb,strand,nb_gap]#,a_length_gap]
                for gen_data in l_gen :
                    if len(gen_data)!=10:
                        print(len(gen_data),'\n')
                        print(gen_data)
                    tab.loc[len(tab.index)] = gen_data
    
    #only keep LCBs longer than 50bp 
    tab=tab.loc[tab['lcb_length'] >= 50]
    
    #sort our df by gen_index then by starting index of nt of each lcb
    l_sort_index = ['gen_index','start_lcb']
    tab = tab.sort_values(l_sort_index,ascending=[True]*len(l_sort_index))
    
    #add a column with the relative nb of lcb
    l_lcb_self_index = list(range(1,nb_lcb+1)) * 2
    
    tab.insert(len(list(tab.columns)), 'lcb_self_index', l_lcb_self_index)
    
    return sv(tab)

#adding sv tag to every sv
def sv(tab):
    
    gen_study = tab.loc[tab['gen_index'] == 2]
    #print(len(gen_study))
    #l_index[0] : index lcb vs reference genome
    #l_index[1] : index lcb for curent genome
    
    #nb_lcb column refer to lcb index of ref genome
    l_nb_lcb = gen_study['nb_lcb'].tolist()
    
    #lcb_index of the genome study
    l_self = gen_study['lcb_self_index'].tolist()
    #print('n_lcb',l_nb_lcb)
    # print('n_self',l_self)
    
    #strand column is only in ref genom so gen_index ==1 and '-' = inversion and '+' = nothing
    l_strand = tab.loc[tab['gen_index'] == 1]['strand'].tolist()

    #inversion, translocation, same
    #set 2 list to store sv at right index for each lcb nb for each genome
    #len(gen_study) = total nb of lcb (same in both genome)
    
    #list of sv for genome study
    l_sv = ['']*len(gen_study)
    #list of sv in ref genome
    l_sv_ref = ['']*len(gen_study)
    
    '''
    to know each sv if there is an insertion of not we are going to check lcb_slef_index vs nb_lcb
    example : 
        nb_lcb          =   1 2 7 3 4 5 6 9 8 
        lcb_self_index  =   1 2 3 4 5 6 7 8 9
        
        index 
        
        nb_lcb[0] (1) = lcb_self_index[0] (1) -> nothing
        
        nb_lcb[1] (2) = lcb_self_index[1] (2) -> nothing 
        nb_lcb[2] (7) = lcb_self_index[1] (3) -> insertion
        but next we will have  nb_lcb[3] (3) = lcb_self_index[3] (4)
        that not the right index for lcb_self_index to compare with nb_lcb because of the insertion
        so we need to adjust the index of lcb_self_index
        to find the insertion in lcb_self_index we are adding a varaible adjust that go up by one if the insertion is from the right to left and retrive one if it is from left to right
        
        if there is an insetion we don't add + 1 to the index only to this when we have an insertion after we keep go up by one for the index
        
        nb_lcb[3] (3) = lcb_self_index[2] (3) -> nothing  
        
        until we found the insertion in lcb_self_index :
            nb_lcb[4] (4) = lcb_self_index[3] (4) -> nothing, lcb_self_index[3] = 4 so 4 + adjust = 4+1 = 5 != 7 
            nb_lcb[5] (5) = lcb_self_index[4] (5) -> nothing    same
            nb_lcb[6] (6) = lcb_self_index[5] (6) -> nothing    lcb_self_index[5] = 6 so 6 + adjust = 6+1 = 7 = 7 we found our sv
            
            we found our sv so
            adjust = 0
            idx += 1 to adjust index that we had stopped going up by one when we found the sv previously
            
            
    '''
    #adjust variable
    adjust = 0
    #set up index_change
    index_change  = np.inf
    
    idx = 0
    for lcb in range(len(l_nb_lcb)):

        if l_self[idx]==l_nb_lcb[lcb]:

            idx+=1
        elif l_self[idx]!=l_nb_lcb[lcb]:
            index_change = l_nb_lcb[lcb]

            if lcb < len(l_nb_lcb)-1:
                ide = lcb
                if l_nb_lcb[lcb+1] > l_nb_lcb[lcb]:
                    adjust -=1
                    
                elif l_nb_lcb[lcb+1] < l_nb_lcb[lcb]: 
                    adjust +=1
                    
                    #when we have 2 sv side by side it is the second one that is the real sv not the 1st one when there are not side by side
                    if l_nb_lcb[lcb+1]+1 == l_nb_lcb[lcb]:
                        ide+=1
            
            l_sv[ide]='INS'
            l_sv_ref[l_nb_lcb[ide]-1]='INS'


        if (l_nb_lcb[lcb] + adjust) == index_change:
            idx +=1
            adjust = 0
            

    #check inversion if there '-' it's an inversion
    for strand in range(len(l_strand)):
        
        if l_strand[strand] == '-':
            index = l_nb_lcb.index(strand+1)
            #print(index,'zoehfzeuiohfzeiou')
            if l_sv[index] == '':
                l_sv[index] = 'INV'
            else:
                l_sv[index] = l_sv[index] + '_INV'
            
            if l_sv_ref[strand] == '':
                l_sv_ref[strand] = 'INV'
            else:
                l_sv_ref[strand] = l_sv_ref[strand] + '_INV'
                
    #create the column sv
    tab['sv']=''
    #add sv to the column
    tab.loc[tab["gen_index"] == 2, "sv"] = l_sv
    tab.loc[tab["gen_index"] == 1, "sv"] = l_sv_ref
    
    return tab

def main():
    
    #check if file is given in input
    if(len(sys.argv)<2):
        print("You have to give the json file with lcb_coordinate given by mauve-parser.js")
        print("USAGE: ", sys.argv[0], '<json_file>')
    else:
        #help
        if sys.argv[1]=='--help' or sys.argv[1]=='--h':
            print("You have to give the json file with lcb_coordinate given by mauve-parser.js")
            print("USAGE: ", sys.argv[0], '<json_file>')
        #check is input file is in json format    
        elif os.path.splitext(sys.argv[1])[1]!='.json':
            print("Please give a json file in input")
        
        else:
            #file_path
            json_file = sys.argv[1]
            
            #extract data to list of list of dico
            data = read_json(json_file)
            
            #data in df 
            df_data = extract_json(data)

            #df to file
            file_name = json_file.split('.xmfa')[0] + '_summary_sorted_ref.csv'
            df_data.to_csv(file_name)
            
            df_data.drop(df_data[df_data.sv == ''].index, inplace=True)
            file_name_sv = json_file.split('.xmfa')[0] + '_summary_sv.csv'
            df_data.to_csv(file_name_sv)
            
            
    return

if __name__ == '__main__':
    main()
