#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:28:22 2022

@author: kot
"""
import pandas as pd
import re
import os
import statistics
import matplotlib.pyplot as plt
import numpy as np
from operator import truediv
'''
    Argument :
        path : path where all csv files are stored
        
    return :
        List of all file name
'''

def get_files(path):
    
    #get all csv files in given path in a list 
    #for root, dirs, files in os.walk(path):
    #    l_files = files.copy()
    
    #regular expression to select only csv files
    extension_ref = re.compile("^[^.].*_summary_sorted_ref.csv")
    extension_sv = re.compile("^[^.].*_summary_sv.csv")
    #keep only csv files in list
    l_summary_files_ref = list(filter(extension_ref.match, os.listdir(path))) #(a,l_files)
    if len(l_summary_files_ref)==0:
        raise Exception("there is 0 ref_csv file in given directory")
    
    l_summary_files_sv = list(filter(extension_sv.match, os.listdir(path))) #(a,l_files)
    if len(l_summary_files_ref)==0:
        raise Exception("there is 0 summary_sv.csv file in given directory")
    return l_summary_files_ref,l_summary_files_sv

'''
function to get all lcb_length for each genome

    Argument:
        file : *_summary_sorted.csv file made by /users/gev/kot/stage/mauve/json_to_csv.py

    return : 
        dico_gen : 
            keys = genome
            values = list of all lcb_length for the specific genome
'''
def extract_lcb(file,data_type):
    
    df = pd.read_csv(file, index_col=0)

    #dico to store data
    dico_gen = {}
    '''
    key ='gen_index'
    if key in df.keys():
        print(df[key])
    else:
        print("{} not a key of dictionary".format(key))
        print(df.keys())
        print(df)'''
    name = []
    if float(df['gen_index'].max()) == 2 :
        
        #find nb of genome
        nb_genome = int(df['gen_index'].max())

        #iterate for each genome
        for genome in range(1, nb_genome + 1):
            
              
            dico_gen[genome] = df.loc[df['gen_index'] == genome][data_type].values.tolist()
            name.append(get_first_val(df,genome))
        
        go = True
    else:
        go = False
    sw = file.split('_')[-4]
    #sw = file.split('.xmfa.json')[0].split('_')[-1]
    return sw, name, dico_gen, go

def get_first_val(df,val):
    try:
        return df.loc[df['gen_index'] == val, 'genome'].iat[0]
    except IndexError:
        return 'no match'



#file = '/users/gev/kot/stage/mauve/seed_weight_test/N2_2013_CB4856_2021_sw_17_summary_sorted.csv'
#dico = extract_lcb(file)
#print(dico)

'''
function to store all lcb length from all genome in dico

    argument :
        path : path of all csv files
'''
def pipe_data(path,data_type):
    
    #gets all csv_file name
    l_summary_files_ref,l_summary_files_sv = get_files(path)
    
    '''
    dico to store data
        key : seed_weight value
        values :
            dico : key = seed_weight value , 
                values : dico : key :   nb_sv : value : nb of sv with this sw
                                        len_all_sv : value : total length of all sv
                                        
    '''
    
    dico_sv = {}
    
    for sv_file in l_summary_files_sv:
        df = pd.read_csv(path + sv_file)
        if len(df)!=0:
            #get seed_weight value
            s_w = int(sv_file.split('_')[-3])
            
            dico_sv[s_w] = {}
            #only go on the ref genome, value of ref genome = 1
            for sv in df.sv.unique():
            
                df_gen1 = df.loc[df['gen_index'] == 1]
                
                dico_sv[s_w][sv] = len(df_gen1.loc[df_gen1['sv'] == sv])
            
            dico_sv[s_w]['nb_sv']=len(df_gen1)
            dico_sv[s_w]['len_all_sv'] = df_gen1['lcb_length'].sum()
    
    
    '''
    dico to store data
        key : seed_weight value
        values :
            dico : key = genome , values = list of all lcb_length
    '''
    dico_sw = {}

    #ieterate for each csv_file
    for file in l_summary_files_ref:
        
        #gets seed_weight value and dico genome and their lcb length
        sw, l_name, dico_gen,go = extract_lcb(path + file,data_type)
        
        if go :
            if list(dico_gen.values()) != [[0], [0]]:
                #add data in dico_sw
                dico_sw[float(sw)] = dico_gen
                dico_sw = dict(sorted(dico_sw.items()))

    return dico_sw, l_name, dico_sv

def merge_repeat(path,d_type):
    start_sw = 17
    end_sw = 31
    
    d_all={}
    d_box = {}
    d_nb_lcb={}
    
    for data_type in d_type:
        d_all[data_type]={}
        if data_type == 'lcb_length':
            d_sv_all={}
            d_sv_al={}
        for seed_weight in range(start_sw, end_sw+1):
            d_all[data_type][seed_weight]={}
            d_sv_all[seed_weight]=[]
            d_sv_al[seed_weight]=[]
            if data_type=='lcb_length':
                
                d_nb_lcb[int(seed_weight)] = {1:[],2:[]}
            for gen in range(1,3):
                d_all[data_type][seed_weight][gen]=[]
    
    for i in range(1,101):
        d_box[i]={}
        print(i)
        path_data = path + str(i) +'/'
        for data_type in d_type:
            
            dico_sw,name,dico_sv = pipe_data(path_data,data_type)
            d_box[i][data_type]=dico_sw
            for sw in list(dico_sw.keys())[1:]:
                
                if data_type=='lcb_length':
                    d_sv_all[int(sw)].append(dico_sv[int(sw)]['len_all_sv'])
                    
                d_sv_al[int(sw)].append(dico_sv[int(sw)]['nb_sv'])
                for ge in list(dico_sw[sw].keys()):
                    d_all[data_type][int(sw)][ge].extend(dico_sw[sw][ge])
                    if data_type=='lcb_length':
                        d_nb_lcb[int(sw)][ge].append(len(dico_sw[sw][ge]))
    
    return d_all, d_sv_all, d_nb_lcb





def plot_data_stat(l_plot, labels, label_data, data_type, genome, name, d_nb_lcb, path_output,repeat):
    
    for plot in range(len(l_plot)) :
                
        if plot == 0 :
            plt.boxplot(l_plot[0], vert=True, patch_artist=True, labels=labels)
            
        elif plot ==3:
            if repeat:
                plt.boxplot(l_plot[3], vert=True, patch_artist=True, labels=list(d_nb_lcb.keys()))
            else:
                plt.plot(labels, l_plot[plot],'*')
        
        else:
            #print(len(l_plot[plot]))
            #print(l_plot[plot])
            plt.plot(labels, l_plot[plot],'*')
            
        title = label_data[plot]+'_'+ data_type + '\n' + name[genome]
        plt.title(title)
        if plot<=2:
            plt.yscale('log')
        plt.xticks(rotation=35)
        plt.xlabel("Seed weight value")
        plt.ylabel(data_type)
        if repeat :
            plt.savefig(path_output + title + '_repeat.png',bbox_inches="tight")
        else:
            plt.savefig(path_output + title + '.png',bbox_inches="tight")
        plt.show()
        plt.close()
    return 





def plot_data(dico_sw,data_type,name,d_nb_lcb,path_output,repeat):

    #get nb of genome
    #nb_gen = max(list(list(dico_sw.values())[0].keys()))
    nb_gen = 2
    #store list of length lcb
    list_of_all_lcb_len_data = []
    #list of seed weight
    labels = list(dico_sw.keys())
    print(labels)
    #add to l_data list of length_lcb
    for genome in range(nb_gen):

        #extract only list of length lcb for each genome
        list_of_values = list(map(lambda d: d[genome+1], dico_sw.values()))

        list_of_all_lcb_len_data.append(list_of_values)
        
        avg_data = []
        nb_lcb = []
        med_data = []
        
        if repeat:
            list_of_nb_lcb = list(map(lambda d: d[1], d_nb_lcb.values()))
            for n_lcb in list_of_nb_lcb:
                #nb_lcb.append(statistics.median(i))
                nb_lcb.append(n_lcb)
        print(nb_lcb)
        data = list_of_all_lcb_len_data[genome]
        
        for length in list_of_all_lcb_len_data[genome]:
            avg_data.append(statistics.mean(length))
            med_data.append(statistics.median(length))
            if not repeat : 
                nb_lcb.append(len(length))
            
        
        
        #fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
        l_plot = [data,avg_data,med_data,nb_lcb]

        label_data = ["boxplot", "average", "median", "number_of_lcb"]#,"length_all_lcb"]
        
        
        plot_data_stat(l_plot, labels, label_data, data_type, genome, name, d_nb_lcb, path_output,repeat)
        
        
        '''
        for plot in range(len(l_plot)) :
            
            if plot == 0 :
                plt.boxplot(l_plot[0], vert=True, patch_artist=True, labels=labels)
                
            elif plot ==3:
                plt.boxplot(l_plot[3], vert=True, patch_artist=True, labels=list(list_of_nb_lcb.keys()))
            
            else:
                #print(len(l_plot[plot]))
                #print(l_plot[plot])
                plt.plot(labels, l_plot[plot][1:],'*')
                
            title = label_data[plot]+'_'+ data_type + '\n' + name[genome]
            plt.title(title)
            if plot<=2:
                plt.yscale('log')
            plt.xticks(rotation=35)
            plt.xlabel("Seed weight value")
            plt.ylabel(data_type)
            plt.savefig('/users/gev/kot/stage/mauve/seed_weed_ws220_cb4856/plot/' + title + '.png',bbox_inches="tight")
            plt.show()
            plt.close()
            '''
    return

def plot_stat_sv(dico_sw, dico_sv,path_output,repeat):
    
    
    if not repeat :
    
        dict_sv = {}
        l_sv_mauve = ['INS','INV','INS_INV','nb_sv','len_all_sv']
        for sv_mauve in l_sv_mauve:
            dict_sv[sv_mauve] = []
        
        l_sw = list(dico_sv.keys())
        for i in range(len(l_sw)):
            for sv in l_sv_mauve:
                if sv in dico_sv[l_sw[i]].keys():
                    dict_sv[sv].append(dico_sv[l_sw[i]][sv])
                else:
                     dict_sv[sv].append(0)
        '''
        for sv in l_sv_mauve:
            plt.plot(list(dico_sv.keys()),dict_sv[sv],'*')
            title = 'number of ' + sv + ' in N2_WS220'
            plt.title(title)
            plt.xticks(rotation=45)
            plt.xlabel("Seed weight value")
            plt.ylabel('number of ' + sv)
            plt.savefig(path_output + title + '.png',bbox_inches="tight")
            plt.show()
            plt.close()
        '''
        l_d = []
        l_sv = []
        list_of_values = list(map(lambda d: d[1], dico_sw.values()))
        print(len(list_of_values))
        for i in range(len(list_of_values)):
            l_d.append(sum(list_of_values[i]))
            l_sv.append(len(list_of_values[i]))
        print(dict_sv['len_all_sv'])
        print(l_d)
        print(len(dict_sv['len_all_sv']))
        print(len(l_d))
        #plt.plot(list(dico_sv.keys()),list(map(truediv, dict_sv['all_sv'],l_plot[-1])),'*')
        labels = list(dico_sv.keys())
        sv_length_ratio = list(map(truediv, dict_sv['len_all_sv'],l_d))
        #a.pop(5)
        #b.pop(5)
        #a.pop(11)
       # b.pop(11)
        #plt.plot(list(dico_sv.keys()),list(map(truediv, dict_sv['all_sv_total_length'],l_d)),'*')
        plt.plot(labels,sv_length_ratio,'*')
        title = 'percentage of sv in all lcb length'
        plt.title(title)
        plt.xticks(rotation=35)
        plt.xlabel("Seed weight value")
        plt.ylabel('percentage of sv in all lcb_length ' + sv)
        plt.savefig(path_output + title + '.png',bbox_inches="tight")
        plt.show()
        plt.close()

    else:
        l_all_lcb = []
        l_all_sv = []
        l_sv = []
        list_of_values = list(map(lambda d: d[1], dico_sw.values()))
        print(len(list_of_values))
        for i in range(len(list_of_values)):
            l_all_lcb.append(sum(list_of_values[i]))
        for j in dico_sv.keys():
            l_all_sv.append(sum(dico_sv[j]))
            l_sv.append(len(dico_sv[j]))
        print(l_all_sv)
        print(l_all_lcb)
        print(len(l_all_sv))
        print(len(l_all_lcb))
        #plt.plot(list(dico_sv.keys()),list(map(truediv, dict_sv['all_sv'],l_plot[-1])),'*')
        labels = list(dico_sv.keys())
        sv_length_ratio = list(map(truediv, l_all_sv,l_all_lcb))
        #a.pop(5)
        #b.pop(5)
        #a.pop(11)
       # b.pop(11)
        #plt.plot(list(dico_sv.keys()),list(map(truediv, dict_sv['all_sv_total_length'],l_d)),'*')
        plt.plot(labels,sv_length_ratio,'*')
        title = 'percentage of sv in all lcb length'
        plt.title(title)
        plt.xticks(rotation=35)
        plt.xlabel("Seed weight value")
        plt.ylabel('percentage of sv in all lcb_length ')
        plt.savefig(path_output + title + '_repeat.png',bbox_inches="tight")
        plt.show()
        plt.close()

    return

'''
for i in range(1,101):
    
    path_output = '/users/gev/kot/stage/mauve/seed_weed_ws220_cb4856/plot/'
    path = '/users/gev/kot/stage/mauve/repeat/seed_weight_N2_WS220_vs_CB4856_19_repeat_24/'# + str(i) + '/'
    d_type = ['lcb_length','nb_gap','average_length_gap']
    repeat = False
    for data_type in d_type:
        dico_sw,name,dico_sv = pipe_data(path,data_type)
        plot_data(dico_sw,data_type,name,dico_sv,path_output,repeat)
        if data_type =='lcb_length':
            plot_stat_sv(dico_sw, dico_sv,path_output,repeat)

'''



#####################repeat#######################
path_output = '/users/gev/kot/stage/mauve/seed_weed_ws220_cb4856/plot/'

path = '/export/home1/users/gev/kot/stage/mauve/repeat/seed_weight_N2_WS220_vs_CB4856_19_repeat_'
d_type = ['lcb_length','nb_gap','average_length_gap']
#d_sw, d_sv, d_n_lcb = merge_repeat(path,d_type)
repeat=True
name = ['N2_WS220','CB_4856_2019']
for data_type in d_type:
    plot_data(d_sw[data_type],data_type,name,d_n_lcb,path_output,repeat)
    if data_type == 'lcb_length':
        plot_stat_sv(d_sw[data_type], d_sv,path_output,repeat)
      
   
        
'''
for i in d_sv_all.keys():
    d_sv_all[i] = {'all_sv_total_length': sum(d_sv_all[i])}

list_of_values = list(map(lambda d: d[1], d_all['lcb_length'].values()))
l_d = []

for i in range(len(list_of_values)):
    l_d.append(list_of_values[i])
len(l_d[0])
for data_type in d_type:       
    plot_data(d_all[data_type],data_type,name,d_sv_all,a)







print(d_all.keys())
d_all['lcb_length']
d_res={}
for i in range(16,32):
    d_res[i]={}
    for data_type in d_type: 
        d_res[i][data_type]={}
        for g in range(1,3):
            d_res[i][data_type][g]={}
            for t in ['med','nb_lcb']:
                d_res[i][data_type][g][t]=[]



for repeat in d_box.keys():
    
    for d_type in d_box[repeat]:
        #get nb of genome
        nb_gen = max(list(list(d_box[repeat][d_type].values())[0].keys()))
        #store list of length lcb
        l_data = []
        #list of seed weight
        labels = list(d_box[repeat][d_type].keys())
        #add to l_data list of length_lcb
        for genome in range(nb_gen):
        
            #extract only list of length lcb for each genome
            list_of_values = list(map(lambda d: d[genome+1], d_box[repeat][d_type].values()))
        
            l_data.append(list_of_values)
        #print(nb_gen)
        #iterate for each genome 
            avg_data = []
            nb_lcb = []
            med_data = []
            data = l_data[genome]
            
            for length in l_data[genome]:
                avg_data.append(statistics.mean(length))
                med_data.append(statistics.median(length))
                nb_lcb.append(len(length))
            
            for s in range(len(med_data)):
            
                d_res[s+16][d_type][genome+1]['med'].append(med_data[s])
                d_res[s+16][d_type][genome+1]['nb_lcb'].append(nb_lcb[s])

#plt.boxplot(l_plot[0], vert=True, patch_artist=True, labels=list(d_res.keys()))
for gen in range(1,3):
    for d_t in ['lcb_length','nb_gap','average_length_gap']:
        med_data = []
        m_data =  []
        for i in d_res.keys():
            med_data.append(d_res[i][d_t][gen]['med'])
            q3, q1 = np.percentile(d_res[i][d_t][gen]['med'], [75 ,25])
            iqr = q3 - q1
            m_data.append(iqr)
        plt.boxplot(med_data[1:], vert=True, patch_artist=True, labels=list(d_res.keys())[1:])
        name = d_t + '_' +str(gen)
        plt.title(name)
        plt.show()
        plt.close()
        print(m_data[1:])
        plt.plot(list(d_res.keys())[1:],m_data[1:],'*')
        name = d_t + '_' +str(gen) + 'iqr'
        plt.title(name)
        #plt.yscale('log')
        plt.show()
        plt.close()



for gen in range(1,3):
    med_data = []
    m_data =  []
    for i in a.keys():
        med_data.append(a[i][gen])
        q3, q1 = np.percentile(a[i][gen], [75 ,25])
        iqr = q3 - q1
        m_data.append(iqr)
    plt.boxplot(med_data, vert=True, patch_artist=True, labels=list(a.keys()))
    name = d_t + '_' +str(gen)

    plt.title(name)
    plt.show()
    plt.close()
    print(m_data)
    plt.plot(list(a.keys()),m_data,'*')
    name = d_t + '_' +str(gen) + 'iqr'
    plt.title(name)
    #plt.yscale('log')
    plt.show()
    plt.close()



m_data =  []
for i in d_sv_al.keys():
    med_data.append(d_sv_al[i][gen])
    q3, q1 = np.percentile(a[i][gen], [75 ,25])
    iqr = q3 - q1
    m_data.append(iqr)
plt.boxplot(med_data, vert=True, patch_artist=True, labels=list(a.keys()))
name = d_t + '_' +str(gen)

plt.title(name)
plt.show()
plt.close()
print(m_data)
plt.plot(list(a.keys()),m_data,'*')
name = d_t + '_' +str(gen) + 'iqr'
plt.title(name)
#plt.yscale('log')
plt.show()
plt.close()


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
            file_name = json_file.split('.')[0] + '_summary_sorted.csv'
            df_data.to_csv(file_name)
            
    return

if __name__ == '__main__':
    main()
'''