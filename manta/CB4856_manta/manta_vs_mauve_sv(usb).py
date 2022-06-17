#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:34:29 2022

@author: kot
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
def plot_manta_vs_mauve (manta_file, mauve_file):
    
    manta_csv = pd.read_csv(manta_file,sep=',')
    manta_csv = manta_csv.iloc[: , 1:]
    print(len(manta_csv['CHROM']=='CHROMOSOME_I'))
    
    mauve_csv = pd.read_csv(mauve_file,sep=',')
    mauve_csv = mauve_csv.iloc[: , 1:]
    print(mauve_csv)
    
    l_chr_manta = manta_csv.CHROM.unique()
    l_chr_mauve = mauve_csv.CHROM.unique()
    
    #l_soft = ['manta','mauve']
    #data_soft = ['dico_chrom','l_chrom_soft','l_svtype','csv_file']
    dico= {}

    dico["manta"]=[{},l_chr_manta, manta_csv]
    dico["mauve"]=[{},l_chr_mauve, mauve_csv]
    
    for software in dico.keys():
        for chrom in dico[software][1]:
            dico[software][0][chrom]={}
            
            pd_chr = dico[software][2].loc[dico[software][2]['CHROM'] == chrom]
            l_sv = pd_chr.SVTYPE.unique()
            
            for sv in l_sv:
                dico[software][0][chrom][sv]= pd_chr.loc[pd_chr['SVTYPE'] == sv][['POS','END']].to_numpy()
    print(dico["mauve"])
    return dico
  
'''          
def dico_chrom(f_chrom_size):
    l_pos = []
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
    with open(f_chrom_size,'r') as f:
        
        for line in f :
            l = line.split(sep='\t')
            l_pos.append(int(l[1]))
    pos_chrom = np.cumsum(l_pos)
    
    dico = {}
    for i in range(len(l_pos)):
        dico[chrom[i]] = pos_chrom[i]
    
    return dico
'''




def plot(dico):
    
    '''
    MAUVE
    INV = b
    INS = r
    INS_INV = darkviolet
    '''
    mauve_sv=['INV','INS','INS_INV']
    
    mauve_colo = ['b','r','darkviolet']
    
    '''
    MANTA
    INV = dodgerblue
    INS = lightcoral
    DEL = g
    BND = black
    DUP = fuchsia
    '''
    
    manta_sv = ['INV','INS','DEL','BND','DUP']
    #soft_sv = [mauve_sv,manta_sv]
    
    
    manta_colo = ['dodgerblue','lightcoral','g','black','fuchsia']
    #soft_colo =[mauve_colo,manta_colo]
    
    dico_col = {}
    dico_col["manta"]= [{},manta_sv,manta_colo]
    dico_col["mauve"]=[{},mauve_sv,mauve_colo]
    for soft in dico_col.keys():
        l_sv = dico_col[soft][1]
        for sv in range(len(l_sv)):
            dico_col[soft][0][l_sv[sv]] = dico_col[soft][2][sv]
    
    #print(dico_col)
    '''
    colo=[]
    n = 10
    
    for i in range(n):
        colo.append('#%06X' % randint(0, 0xFFFFFF))
    '''    
    
    
    nb_sv = 0
    l=[]
    for software in dico.keys():
        nb_sv += len(dico[software][2])
        l.append(dico[software][1])
    
    li = l[0].tolist() + list(set(l[1].tolist()) - set(l[0].tolist()))

    l_nb_sv = np.linspace(0,1,num=nb_sv)
    np.random.shuffle(l_nb_sv)
    n = 0
    
    chrom = ['CHROMOSOME_I','CHROMOSOME_II','CHROMOSOME_III','CHROMOSOME_IV','CHROMOSOME_V','CHROMOSOME_X','CHROMOSOME_MtDNA']
    d_nb_sv_manta ={}
    d_nb_sv_mauve ={}
    l_length_chr = [15072423,15279345,13783700,17493793,20924149,17718866,13794]

    
    for i in range(len(chrom)):
        d_nb_sv_manta[chrom[i]]=[0]*(l_length_chr[i])
        d_nb_sv_mauve[chrom[i]]=[0]*(l_length_chr[i])
        #for j in range(0,l_length_chr[i],10000):
            #d_nb_sv_manta[chrom[i]][j]=0
            #d_nb_sv_mauve[chrom[i]][j]=0
    
    
    
    #fig, ax = plt.subplots()
    l_label =[[],[],[]]
    for chrom_name in li:
       # nb_chr = 0
        
        na=0
        fig, ax = plt.subplots()
        
        for software in dico.keys():
            if chrom_name in dico[software][0].keys():
            #chrom_name = dico[software][1][chrom]
                for svtype in range(len(list(dico[software][0][chrom_name].keys()))):
                    if software=='mauve':
                        print(list(dico[software][0][chrom_name].keys()))
                    col_sv = dico_col[software][0][list(dico[software][0][chrom_name].keys())[svtype]]
                    #print(col_sv)
                    for coord_sv in dico[software][0][chrom_name][list(dico[software][0][chrom_name].keys())[svtype]]:
                        na+=1
                        if coord_sv[1]!=-1: #color=colo[svtype]
                           if software == 'manta':
                               d_nb_sv_manta[chrom_name][coord_sv[0]]+=1
                               ax.plot(coord_sv, [l_nb_sv[n]]*2, color=col_sv, linestyle='-', marker='', linewidth=2.0)
                               leng = (coord_sv [1] - coord_sv[0])
                               if leng>3000000:
                                   print('pos : ', coord_sv)
                           else:
                                d_nb_sv_mauve[chrom_name][coord_sv[0]]+=1
                                
                                ax.plot(coord_sv, [l_nb_sv[n]]*2, color=col_sv, linestyle='-', marker='', linewidth=2.0)
                                
                                #ax.plot(coord_sv, [l_nb_sv[n]]*2, color='red', linestyle='-', marker='', linewidth=3.0)
                                #for i in range(coord_sv[0],coord_sv[1]+1,siz):
                                    #d_nb_sv_mauve[chrom_name][i]+=1
                            #l_label[1].append((coord_sv[1]-coord_sv[0])/2 + coord_sv[0])
                            #l_label[1].append(coord_sv)
                            
                        else:
                            ax.plot(coord_sv[0], l_nb_sv[n], '*',color=col_sv)
                            #l_label[1].append(coord_sv[0])
                            d_nb_sv_manta[chrom_name][coord_sv[0]]+=1
                        
                        
                        
                        l_label[0].append(software + '_' + list(dico[software][0][chrom_name].keys())[svtype])
                        l_label[2].append(l_nb_sv[n])
                        
                        n+=1
        
        '''
        tag = [[],[],[]]
        wide = []
        size = 1000
        for i in range(len(l_label[0])):
            if l_label[0][i].split(sep='_')[0]=='mauve':
                wide.append(l_label[1][i])
        for j in range(len(l_label[1])):
            k=0
            while k < len(wide) or l_label[1][j] + size :
                if l_label[1][j] + size 
                
        
        '''
        '''
        ###############LABEL ####################################                
        for label, x, y in zip(l_label[0], l_label[1], l_label[2]):
            if label.split(sep='_')[0]=='mauve':
                c = 'green'
            else:
                c= 'yellow'
            plt.annotate(
                label,
                xy=(x, y), xytext=(-20, 20),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc=c, alpha=0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        '''
        
        
        pop_a = mpatches.Patch(color='b', label='Mauve_INV')
        pop_b = mpatches.Patch(color='r', label='Mauve_INS')
        pop_c = mpatches.Patch(color='darkviolet', label='Mauve_INS_INV')
        pop_d = mpatches.Patch(color='dodgerblue', label='Manta_INV')
        pop_e = mpatches.Patch(color='lightcoral', label='Manta_INS')
        pop_f = mpatches.Patch(color='g', label='Manta_DEL')
        pop_g = mpatches.Patch(color='black', label='Manta_BND')
        pop_h = mpatches.Patch(color='fuchsia', label='Manta_DUP')
        
        plt.legend(handles=[pop_a,pop_b,pop_c,pop_d,pop_e,pop_f,pop_g,pop_h],bbox_to_anchor=(1,1), loc="upper left")
        
        print('nb_sv ', chrom_name , ' : ', na)
        #plt.xlim([0, l_length_chr[chrom]])
        plt.title(chrom_name)
        #plt.legend()
        plt.show()
        #plt.close()
        
        '''
        #x =  d_nb_sv_manta[chrom_name]

        #x_nb_manta = list(range(1,len(d_nb_sv_manta[chrom_name])+1))
        #x_nb_mauve = list(range(1,len(d_nb_sv_manta[chrom_name])+1))
        
        
        
        y_nb_manta = np.cumsum(d_nb_sv_manta[chrom_name])
        #y_nb_manta = [i for i in y_nb_manta if i != 0]

        y_nb_mauve = np.cumsum(d_nb_sv_mauve[chrom_name])
        #y_nb_mauve = [i for i in y_nb_mauve if i != 0]
        x_nb_manta=[]
        x_nb_mauve=[]
        for i in range(len(y_nb_manta)) : 
            if y_nb_manta[i]!=0:
                x_nb_manta.append(i+1)
        for i in range(len(y_nb_mauve)) : 
            if y_nb_mauve[i]!=0:
                x_nb_mauve.append(i+1)
                
        #y_nb_manta = [i for i in y_nb_manta if i != 0]        
        #y_nb_mauve = [i for i in y_nb_mauve if i != 0]
        plt.plot(y_nb_manta)#,x_nb_manta,'*')
        plt.plot(y_nb_mauve)#,x_nb_mauve,'+')
        
        
        
        plt.title(chrom_name)
        plt.show()
        
        #nb_chr =+ 1
        '''
        
    return

#manta_file = '/media/julien/8BC0-9F2B/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv_compare_sv_qual_0_hom_False.csv'
#manta_file = '/media/julien/8BC0-9F2B/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv_compare_sv_qual_900_hom_False.csv'
#manta_file = '/media/julien/8BC0-9F2B/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv_compare_sv_qual_0_hom_True.csv'
manta_file = '/media/kot/8BC0-9F2B/manta/CB4856_manta/CB4856.A6330/results/variants/diploidSV_inv_compare_sv_qual_900_hom_True.csv'
#mauve_file = '/users/gev/kot/stage/manta/CB4856_manta/sv_cb4856_mauve.csv'
mauve_file = '/media/kot/8BC0-9F2B/Mauve/seed_ws220/csv/csv/WS220_CB4856_2019_sw_20_summary_sv_compare_sv.csv'
#print(plot_manta_vs_mauve(manta_file, mauve_file)['mauve'])
plot(plot_manta_vs_mauve(manta_file, mauve_file))


'''
import plotly.graph_objects as px
import plotly.express as go
import numpy as np
  
df = go.data.tips()
  
x = df['total_bill']
y = df['tip']
  
plot = px.Figure(data=[px.Scatter(
    x=x,
    y=y,
    mode='markers',)
])
  
plot.update_layout(
    xaxis=dict(
        rangeselector=dict(
            buttons=list([
                dict(count=1,
                    step="day",
                    stepmode="backward"),
            ])
        ),
        rangeslider=dict(
            visible=True
        ),
    )
)
  
plot.show()


'''




'''    for software in dico.keys():
        for chrom in range(len(dico[software][1])):
            fig, ax = plt.subplots()
            chrom_name = dico[software][1][chrom]
            for svtype in range(len(list(dico[software][0][chrom_name].keys()))):
                for sv in dico[software][0][chrom_name][list(dico[software][0][chrom_name].keys())[svtype]]:
                    if sv[1]!=-1:    
                        ax.plot(sv, [l_nb_sv[n]]*2, color=colo[svtype], linestyle='-', marker='')
                    else:
                        ax.plot(sv[0], l_nb_sv[n], color=colo[svtype], linestyle='-', marker='',linewidth=5)
                    n+=1
            #plt.xlim([0, l_length_chr[chrom]])
            plt.title(chrom_name)
            plt.show()
            plt.close() '''
            
'''

manta_csv = pd.read_csv(manta_file,sep=',')
manta_csv = manta_csv.iloc[: , 1:]
print(manta_csv.loc[manta_csv['CHROM'] == 'CHROMOSOME_I'])


manta_csv.loc[((manta_csv['CHROM'] == 'CHROMOSOME_I') & (manta_csv['SVTYPE'] == 'INV'))]['POS'].tolist()


mauve_csv.plot('POS','END')

x1 = [1, 100]
x2 = [1, 10**4]    
ys = [1, 1]

fig, ax = plt.subplots()
ax.plot(x1, ys, 'gray', linestyle=':', marker='')
ax.plot(x2, ys, 'black', linestyle='--', marker='')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([1, 10, 100, 1000, 10000])
ax.set_yticks([10**-4, 10**-3, 10**-2, 10**-1, 10**0, 10**1])
ax.get_yaxis().get_major_formatter().labelOnlyBase = False
ax.get_xaxis().get_major_formatter().labelOnlyBase = False
plt.show()





import random
  
# list of items
List = [10, 20, 30, 40, 50, 40,
        30, 20, 10]
  
# using the sample() method
random.shuffle(List)
  
# displaying random selections from 
# the list without repetition
print(List)
'''