###############################################
##Dmitry Sutormin, 2022##
##Ago-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over intergenic regions (IGRs). Plots this information.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
from matplotlib_venn import venn2, venn3, venn3_circles
from matplotlib import cm as cm
import collections
from collections import OrderedDict
import pandas as pd
from pandas import DataFrame


#Path to the working directory.
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Metaintergenic_analysis\IGR_polyshed\\'

#Half-window width will be used to smooth signal.
Sm_window=20
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_IGR={'Lr_5.5_cip_F ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_F +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_F -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_F --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_R ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_R +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_R -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_cip_R --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_cip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_F ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_F +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_F -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_F --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_R ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_R +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_R -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_cip_R --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_cip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_F ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_F +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_F -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_F --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_R ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_R +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_R -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_cip_R --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_cip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_F ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_F +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_F -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_F --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_R ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_R +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_R -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_cip_R --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_cip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_F ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_F +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_F -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_F --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_R ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_R +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_R -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_5.5_nocip_R --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_5.5_nocip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_F ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_F +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_F -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_F --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_R ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_R +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_R -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Lr_12.5_nocip_R --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Lr_12.5_nocip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_F ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_F +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_F -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_F --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_R ++'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_R +-'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_R -+'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_5.5_nocip_R --'     :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_5.5_nocip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_F ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_F_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_F +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_F_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_F -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_F_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_F --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_F_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_R ++'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_R_over_++IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_R +-'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_R_over_+-IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_R -+'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_R_over_-+IGR_50_1000_all_width_200bp_gb_100bp.wig',
                      'Se_12.5_nocip_R --'    :    PWD + 'Signal_of_TUs_wig\IGR_50_1000_all\Signal_Se_12.5_nocip_R_over_--IGR_50_1000_all_width_200bp_gb_100bp.wig', 
                      }

#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="IGR"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Figures\Plot_combinations\\'
Dir_check_create(Out_path)


#######
#Parses WIG file with FE over TUs.
#######

def wig_FE_over_genes_parsing(name, wigfile):
    print('Now is processing: ' + str(name) + ' ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] in ['track']:
            ww_l=line[2].split('=')[1].rstrip('"').lstrip('"').split('_')
            win_width=int(ww_l[0])
            length=int(ww_l[1])
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(f'Window width: {win_width}, length of TU: {length}')
    return NE_values, win_width, length


#######
#Returns smoothed tracks.
#######

def Smoothing(ends, window):
    smoothed=[]
    #Calculating the value for the first position
    sm=0.0
    window_float=float(window)
    sm+=np.mean(ends[:2*window])
    smoothed.append(sm)
    #Calculating values for the part of the array remains
    for i in range(len(ends)-2*window):
        sm+=(ends[i+2*window]-ends[i])/(2*window_float)
        smoothed.append(sm)
    return smoothed


#######
#Plot the signal for different groups of genes together.
#######

def plot_intergenic_FE_TUs_groups(wig_in_dict, sm_window, output_path, set_name, set_type):
    
    #Number of intervals in sets.
    IGR_sets_size={'Lr_5.5_cip_F ++'  : 690,
                   'Lr_5.5_cip_F +-'  : 325,
                   'Lr_5.5_cip_F -+'  : 603,
                   'Lr_5.5_cip_F --'  : 735,
                   'Lr_5.5_cip_R ++'  : 690,
                   'Lr_5.5_cip_R +-'  : 325,
                   'Lr_5.5_cip_R -+'  : 603,
                   'Lr_5.5_cip_R --'  : 735,
                   'Lr_12.5_cip_F ++' : 690,
                   'Lr_12.5_cip_F +-' : 325,
                   'Lr_12.5_cip_F -+' : 603,
                   'Lr_12.5_cip_F --' : 735,
                   'Lr_12.5_cip_R ++' : 690,
                   'Lr_12.5_cip_R +-' : 325,
                   'Lr_12.5_cip_R -+' : 603,
                   'Lr_12.5_cip_R --' : 735,
                   'Se_5.5_cip_F ++'  : 690,
                   'Se_5.5_cip_F +-'  : 325,
                   'Se_5.5_cip_F -+'  : 603,
                   'Se_5.5_cip_F --'  : 735,
                   'Se_5.5_cip_R ++'  : 690,
                   'Se_5.5_cip_R +-'  : 325,
                   'Se_5.5_cip_R -+'  : 603,
                   'Se_5.5_cip_R --'  : 735,
                   'Se_12.5_cip_F ++' : 690,
                   'Se_12.5_cip_F +-' : 325,
                   'Se_12.5_cip_F -+' : 603,
                   'Se_12.5_cip_F --' : 735,
                   'Se_12.5_cip_R ++' : 690,
                   'Se_12.5_cip_R +-' : 325,
                   'Se_12.5_cip_R -+' : 603,
                   'Se_12.5_cip_R --' : 735,
                   'Lr_5.5_nocip_F ++'  : 690,
                   'Lr_5.5_nocip_F +-'  : 325,
                   'Lr_5.5_nocip_F -+'  : 603,
                   'Lr_5.5_nocip_F --'  : 735,
                   'Lr_5.5_nocip_R ++'  : 690,
                   'Lr_5.5_nocip_R +-'  : 325,
                   'Lr_5.5_nocip_R -+'  : 603,
                   'Lr_5.5_nocip_R --'  : 735,
                   'Lr_12.5_nocip_F ++' : 690,
                   'Lr_12.5_nocip_F +-' : 325,
                   'Lr_12.5_nocip_F -+' : 603,
                   'Lr_12.5_nocip_F --' : 735,
                   'Lr_12.5_nocip_R ++' : 690,
                   'Lr_12.5_nocip_R +-' : 325,
                   'Lr_12.5_nocip_R -+' : 603,
                   'Lr_12.5_nocip_R --' : 735,
                   'Se_5.5_nocip_F ++'  : 690,
                   'Se_5.5_nocip_F +-'  : 325,
                   'Se_5.5_nocip_F -+'  : 603,
                   'Se_5.5_nocip_F --'  : 735,
                   'Se_5.5_nocip_R ++'  : 690,
                   'Se_5.5_nocip_R +-'  : 325,
                   'Se_5.5_nocip_R -+'  : 603,
                   'Se_5.5_nocip_R --'  : 735,
                   'Se_12.5_nocip_F ++' : 690,
                   'Se_12.5_nocip_F +-' : 325,
                   'Se_12.5_nocip_F -+' : 603,
                   'Se_12.5_nocip_F --' : 735,
                   'Se_12.5_nocip_R ++' : 690,
                   'Se_12.5_nocip_R +-' : 325,
                   'Se_12.5_nocip_R -+' : 603,
                   'Se_12.5_nocip_R --' : 735,              
                   }    
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=5000
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(positions[0], positions[-1])
       
    
    #Plot FE over genes.
    plt.figure(figsize=(4, 4), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_name=="Lr_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.        
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_R ++']) + np.array(dict_of_wigs['Lr_5.5_cip_F ++']))/(np.array(dict_of_wigs['Lr_5.5_nocip_R ++']) + np.array(dict_of_wigs['Lr_5.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Lr_5.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_R +-']) + np.array(dict_of_wigs['Lr_5.5_cip_F +-']))/(np.array(dict_of_wigs['Lr_5.5_nocip_R +-']) + np.array(dict_of_wigs['Lr_5.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Lr_5.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_R -+']) + np.array(dict_of_wigs['Lr_5.5_cip_F -+']))/(np.array(dict_of_wigs['Lr_5.5_nocip_R -+']) + np.array(dict_of_wigs['Lr_5.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Lr_5.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_R --']) + np.array(dict_of_wigs['Lr_5.5_cip_F --']))/(np.array(dict_of_wigs['Lr_5.5_nocip_R --']) + np.array(dict_of_wigs['Lr_5.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Lr_5.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5":    
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_R ++']) + np.array(dict_of_wigs['Se_5.5_cip_F ++']))/(np.array(dict_of_wigs['Se_5.5_nocip_R ++']) + np.array(dict_of_wigs['Se_5.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Se_5.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_R +-']) + np.array(dict_of_wigs['Se_5.5_cip_F +-']))/(np.array(dict_of_wigs['Se_5.5_nocip_R +-']) + np.array(dict_of_wigs['Se_5.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Se_5.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_R -+']) + np.array(dict_of_wigs['Se_5.5_cip_F -+']))/(np.array(dict_of_wigs['Se_5.5_nocip_R -+']) + np.array(dict_of_wigs['Se_5.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Se_5.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_R --']) + np.array(dict_of_wigs['Se_5.5_cip_F --']))/(np.array(dict_of_wigs['Se_5.5_nocip_R --']) + np.array(dict_of_wigs['Se_5.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Se_5.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5":     
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_R ++']) + np.array(dict_of_wigs['Se_12.5_cip_F ++']))/(np.array(dict_of_wigs['Se_12.5_nocip_R ++']) + np.array(dict_of_wigs['Se_12.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Se_12.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_R +-']) + np.array(dict_of_wigs['Se_12.5_cip_F +-']))/(np.array(dict_of_wigs['Se_12.5_nocip_R +-']) + np.array(dict_of_wigs['Se_12.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Se_12.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_R -+']) + np.array(dict_of_wigs['Se_12.5_cip_F -+']))/(np.array(dict_of_wigs['Se_12.5_nocip_R -+']) + np.array(dict_of_wigs['Se_12.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Se_12.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_R --']) + np.array(dict_of_wigs['Se_12.5_cip_F --']))/(np.array(dict_of_wigs['Se_12.5_nocip_R --']) + np.array(dict_of_wigs['Se_12.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Se_12.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5":     
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_R ++']) + np.array(dict_of_wigs['Lr_12.5_cip_F ++']))/(np.array(dict_of_wigs['Lr_12.5_nocip_R ++']) + np.array(dict_of_wigs['Lr_12.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Lr_12.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_R +-']) + np.array(dict_of_wigs['Lr_12.5_cip_F +-']))/(np.array(dict_of_wigs['Lr_12.5_nocip_R +-']) + np.array(dict_of_wigs['Lr_12.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Lr_12.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_R -+']) + np.array(dict_of_wigs['Lr_12.5_cip_F -+']))/(np.array(dict_of_wigs['Lr_12.5_nocip_R -+']) + np.array(dict_of_wigs['Lr_12.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Lr_12.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_R --']) + np.array(dict_of_wigs['Lr_12.5_cip_F --']))/(np.array(dict_of_wigs['Lr_12.5_nocip_R --']) + np.array(dict_of_wigs['Lr_12.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Lr_12.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
   
                    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS/TE'
    ticks_lables[ticks.index(length)]='TE/TS'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE/TS')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_cip_div_nocip_FR_{win_width}bp_with_IGR_{length}_bp.png', dpi=400, figsize=(4, 4), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(4, 4), dpi=100)
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_name=="Lr_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.        
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_R ++']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_F ++']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_R ++']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Lr_5.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_R +-']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_F +-']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_R +-']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Lr_5.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_R -+']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_F -+']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_R -+']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Lr_5.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_R --']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_F --']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_R --']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Lr_5.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5":    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_R ++']) + np.array(dict_of_wigs_sm['Se_5.5_cip_F ++']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_R ++']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Se_5.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_R +-']) + np.array(dict_of_wigs_sm['Se_5.5_cip_F +-']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_R +-']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Se_5.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_R -+']) + np.array(dict_of_wigs_sm['Se_5.5_cip_F -+']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_R -+']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Se_5.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_R --']) + np.array(dict_of_wigs_sm['Se_5.5_cip_F --']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_R --']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Se_5.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5":     
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_R ++']) + np.array(dict_of_wigs_sm['Se_12.5_cip_F ++']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_R ++']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Se_12.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_R +-']) + np.array(dict_of_wigs_sm['Se_12.5_cip_F +-']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_R +-']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Se_12.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_R -+']) + np.array(dict_of_wigs_sm['Se_12.5_cip_F -+']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_R -+']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Se_12.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_R --']) + np.array(dict_of_wigs_sm['Se_12.5_cip_F --']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_R --']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Se_12.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5":     
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_R ++']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_F ++']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_R ++']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_F ++'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'-> -> ({IGR_sets_size["Lr_12.5_cip_R ++"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_R +-']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_F +-']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_R +-']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_F +-'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'-> <- ({IGR_sets_size["Lr_12.5_cip_R +-"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_R -+']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_F -+']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_R -+']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_F -+'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'<- -> ({IGR_sets_size["Lr_12.5_cip_R -+"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_R --']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_F --']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_R --']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_F --'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'<- <- ({IGR_sets_size["Lr_12.5_cip_R --"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
          
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS/TE'
    ticks_lables[ticks.index(length)]='TE/TS'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE/TS')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    #plot1.axhline(1, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)
    #plot1.set_title(f'{set_name} signal over {set_type}', size=20)   
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_cip_div_nocip_FR_smoothed_{win_width}bp_with_IGR_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(4, 4), transparent=True)   
    plt.show()
    plt.close()    
    return

#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_5.5'
plot_intergenic_FE_TUs_groups(Wig_data_in_dict_IGR, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5'
plot_intergenic_FE_TUs_groups(Wig_data_in_dict_IGR, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5'
plot_intergenic_FE_TUs_groups(Wig_data_in_dict_IGR, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5'
plot_intergenic_FE_TUs_groups(Wig_data_in_dict_IGR, Sm_window, Out_path, Signal_name, Set_type)