###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#signal over symmetric anchors (e.g., GCSs). Plots this information.

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
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Anchor_analysis_polyshed\\'

#Half-window width will be used to smooth signal.
Sm_window=1000
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_GSCs={'Lr_5.5_cip_F_0.01_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_5.5_cip_R_0.01_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_5.5_cip_F_0.05_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_5.5_cip_R_0.05_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_5.5_cip_F_0.01_top'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_5.5_cip_R_0.01_top'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_5.5_cip_F_0.01_bot'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Lr_5.5_cip_R_0.01_bot'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                          
                       'Lr_5.5_nocip_F_0.01_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_5.5_nocip_R_0.01_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_5.5_nocip_F_0.05_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_5.5_nocip_R_0.05_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_5.5_nocip_F_0.01_top'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_5.5_nocip_R_0.01_top'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_5.5_nocip_F_0.01_bot'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Lr_5.5_nocip_R_0.01_bot'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                       
                       'Lr_12.5_cip_F_0.01_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_12.5_cip_R_0.01_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_12.5_cip_F_0.05_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_12.5_cip_R_0.05_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_12.5_cip_F_0.01_top'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_12.5_cip_R_0.01_top'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_12.5_cip_F_0.01_bot'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Lr_12.5_cip_R_0.01_bot'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                          
                       'Lr_12.5_nocip_F_0.01_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_12.5_nocip_R_0.01_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Lr_12.5_nocip_F_0.05_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_12.5_nocip_R_0.05_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Lr_12.5_nocip_F_0.01_top'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_12.5_nocip_R_0.01_top'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Lr_12.5_nocip_F_0.01_bot'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Lr_12.5_nocip_R_0.01_bot'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                       
                       'Se_5.5_cip_F_0.01_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Se_5.5_cip_R_0.01_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Se_5.5_cip_F_0.05_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Se_5.5_cip_R_0.05_all'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Se_5.5_cip_F_0.01_top'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_5.5_cip_R_0.01_top'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_5.5_cip_F_0.01_bot'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Se_5.5_cip_R_0.01_bot'     :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                                                                                                                                 
                       'Se_5.5_nocip_F_0.01_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Se_5.5_nocip_R_0.01_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Se_5.5_nocip_F_0.05_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Se_5.5_nocip_R_0.05_all'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Se_5.5_nocip_F_0.01_top'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_5.5_nocip_R_0.01_top'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_5.5_nocip_F_0.01_bot'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Se_5.5_nocip_R_0.01_bot'   :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                                                                                                                                 
                       'Se_12.5_cip_F_0.01_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Se_12.5_cip_R_0.01_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Se_12.5_cip_F_0.05_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Se_12.5_cip_R_0.05_all'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Se_12.5_cip_F_0.01_top'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_12.5_cip_R_0.01_top'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_12.5_cip_F_0.01_bot'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Se_12.5_cip_R_0.01_bot'    :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_bot_width_15000bp.wig',
                                                                                                                                 
                       'Se_12.5_nocip_F_0.01_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_width_15000bp.wig',
                       'Se_12.5_nocip_R_0.01_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_width_15000bp.wig',
                       'Se_12.5_nocip_F_0.05_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.05_width_15000bp.wig',
                       'Se_12.5_nocip_R_0.05_all'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.05_width_15000bp.wig',
                       'Se_12.5_nocip_F_0.01_top'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_12.5_nocip_R_0.01_top'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_top_width_15000bp.wig',
                       'Se_12.5_nocip_F_0.01_bot'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_bot_width_15000bp.wig',
                       'Se_12.5_nocip_R_0.01_bot'  :  PWD + 'GCSs_as_anchors_15000_width_200_smoothing\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_bot_width_15000bp.wig',                                                                        
                       }

#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="GCSs"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\GCSs_as_anchors_15000_width_200_smoothing\Plot_combinations\\'
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
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(f'Window width: {win_width}')
    return NE_values, win_width


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
#Plot the signal for different groups of samples together.
#Fold enrichment or other ratios.
#######

def plot_FE_over_anchors(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    GSCs_sets={'Lr_5.5_cip_F_0.01_all'     :  1833,
               'Lr_5.5_cip_R_0.01_all'     :  1833,
               'Lr_5.5_cip_F_0.05_all'     :  3828,
               'Lr_5.5_cip_R_0.05_all'     :  3828,
               'Lr_5.5_cip_F_0.01_top'     :  184,
               'Lr_5.5_cip_R_0.01_top'     :  184,
               'Lr_5.5_cip_F_0.01_bot'     :  193,
               'Lr_5.5_cip_R_0.01_bot'     :  193,
                                              
               'Lr_5.5_nocip_F_0.01_all'   :  1833,
               'Lr_5.5_nocip_R_0.01_all'   :  1833,
               'Lr_5.5_nocip_F_0.05_all'   :  3828,
               'Lr_5.5_nocip_R_0.05_all'   :  3828,
               'Lr_5.5_nocip_F_0.01_top'   :  184,
               'Lr_5.5_nocip_R_0.01_top'   :  184,
               'Lr_5.5_nocip_F_0.01_bot'   :  193,
               'Lr_5.5_nocip_R_0.01_bot'   :  193,
                                              
               'Lr_12.5_cip_F_0.01_all'    :  1833,
               'Lr_12.5_cip_R_0.01_all'    :  1833,
               'Lr_12.5_cip_F_0.05_all'    :  3828,
               'Lr_12.5_cip_R_0.05_all'    :  3828,
               'Lr_12.5_cip_F_0.01_top'    :  184,
               'Lr_12.5_cip_R_0.01_top'    :  184,
               'Lr_12.5_cip_F_0.01_bot'    :  193,
               'Lr_12.5_cip_R_0.01_bot'    :  193,
                                              
               'Lr_12.5_nocip_F_0.01_all'  :  1833,
               'Lr_12.5_nocip_R_0.01_all'  :  1833,
               'Lr_12.5_nocip_F_0.05_all'  :  3828,
               'Lr_12.5_nocip_R_0.05_all'  :  3828,
               'Lr_12.5_nocip_F_0.01_top'  :  184,
               'Lr_12.5_nocip_R_0.01_top'  :  184,
               'Lr_12.5_nocip_F_0.01_bot'  :  193,
               'Lr_12.5_nocip_R_0.01_bot'  :  193,
                                              
               'Se_5.5_cip_F_0.01_all'     :  1833,
               'Se_5.5_cip_R_0.01_all'     :  1833,
               'Se_5.5_cip_F_0.05_all'     :  3828,
               'Se_5.5_cip_R_0.05_all'     :  3828,
               'Se_5.5_cip_F_0.01_top'     :  184,
               'Se_5.5_cip_R_0.01_top'     :  184,
               'Se_5.5_cip_F_0.01_bot'     :  193,
               'Se_5.5_cip_R_0.01_bot'     :  193,
                                              
               'Se_5.5_nocip_F_0.01_all'   :  1833,
               'Se_5.5_nocip_R_0.01_all'   :  1833,
               'Se_5.5_nocip_F_0.05_all'   :  3828,
               'Se_5.5_nocip_R_0.05_all'   :  3828,
               'Se_5.5_nocip_F_0.01_top'   :  184,
               'Se_5.5_nocip_R_0.01_top'   :  184,
               'Se_5.5_nocip_F_0.01_bot'   :  193,
               'Se_5.5_nocip_R_0.01_bot'   :  193,
                                              
               'Se_12.5_cip_F_0.01_all'    :  1833,
               'Se_12.5_cip_R_0.01_all'    :  1833,
               'Se_12.5_cip_F_0.05_all'    :  3828,
               'Se_12.5_cip_R_0.05_all'    :  3828,
               'Se_12.5_cip_F_0.01_top'    :  184,
               'Se_12.5_cip_R_0.01_top'    :  184,
               'Se_12.5_cip_F_0.01_bot'    :  193,
               'Se_12.5_cip_R_0.01_bot'    :  193,
                                              
               'Se_12.5_nocip_F_0.01_all'  :  1833,
               'Se_12.5_nocip_R_0.01_all'  :  1833,
               'Se_12.5_nocip_F_0.05_all'  :  3828,
               'Se_12.5_nocip_R_0.05_all'  :  3828,
               'Se_12.5_nocip_F_0.01_top'  :  184,
               'Se_12.5_nocip_R_0.01_top'  :  184,
               'Se_12.5_nocip_F_0.01_bot'  :  193,
               'Se_12.5_nocip_R_0.01_bot'  :  193,                                                                                       
               }    
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=4
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(positions[0], positions[-1])
    
    #Plot FE over genes.
    plt.figure(figsize=(4, 4), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Anchor sets below.
    if set_name=="Lr_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_top']))/(np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.05_all']))/(np.array(dict_of_wigs['Lr_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_top']))/(np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.05_all']))/(np.array(dict_of_wigs['Lr_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_top']))/(np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.05_all']))/(np.array(dict_of_wigs['Se_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_top']))/(np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.05_all']))/(np.array(dict_of_wigs['Se_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
    elif set_name=="Comp_cond":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_all'])) /(np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_all'])),  linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_all'])+ np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_all'])+ np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_all'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_all'])) /(np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_all'])),  linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Se_12.5_cip_F_0.01_all'])+ np.array(dict_of_wigs['Se_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_all'])+ np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        
        
                    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1, 5000).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)   
    #plot1.set_xlim([-1500, 1500])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_cip_div_nocip_FR_{win_width}bp_with_IGR_{length}_bp.png', dpi=400, figsize=(4, 4), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
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
    
    ##Anchor sets below.
    if set_name=="Lr_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_top']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.05_all']))/(np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_top']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.05_all']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_top']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.05_all']))/(np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_all'])), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_bot']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_bot'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_top']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_top'])), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.05_all']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.05_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
    elif set_name=="Comp_cond":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_all'])) /(np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_all'])),  linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_all'])+ np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_all'])+ np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_all'])), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_all'])) /(np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_all'])),  linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_all'])+ np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_all']))/(np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_all'])+ np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_all'])), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
       
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1, 5000).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20) 
    #plot1.set_xlim([-1500, 1500])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_cip_div_nocip_FR_smoothed_{win_width}bp_with_IGR_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(4, 4), transparent=True)   
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_cip_div_nocip_FR_smoothed_{win_width}bp_with_IGR_{length}bp_smoothed_{2*sm_window}bp.svg', dpi=400, figsize=(4, 4), transparent=True)   
    plt.close()  
    
    return

#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_5.5'
plot_FE_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_12.5'
plot_FE_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_5.5'
plot_FE_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_12.5'
plot_FE_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Comp_cond'
plot_FE_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)



#######
#Plot the signal for different groups of anchors together.
#Coverage depth or RPKM.
#######

def plot_RPKM_over_anchors(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    GSCs_sets={'Lr_5.5_cip_F_0.01_all'     :  1833,
               'Lr_5.5_cip_R_0.01_all'     :  1833,
               'Lr_5.5_cip_F_0.05_all'     :  3828,
               'Lr_5.5_cip_R_0.05_all'     :  3828,
               'Lr_5.5_cip_F_0.01_top'     :  184,
               'Lr_5.5_cip_R_0.01_top'     :  184,
               'Lr_5.5_cip_F_0.01_bot'     :  193,
               'Lr_5.5_cip_R_0.01_bot'     :  193,
                                              
               'Lr_5.5_nocip_F_0.01_all'   :  1833,
               'Lr_5.5_nocip_R_0.01_all'   :  1833,
               'Lr_5.5_nocip_F_0.05_all'   :  3828,
               'Lr_5.5_nocip_R_0.05_all'   :  3828,
               'Lr_5.5_nocip_F_0.01_top'   :  184,
               'Lr_5.5_nocip_R_0.01_top'   :  184,
               'Lr_5.5_nocip_F_0.01_bot'   :  193,
               'Lr_5.5_nocip_R_0.01_bot'   :  193,
                                              
               'Lr_12.5_cip_F_0.01_all'    :  1833,
               'Lr_12.5_cip_R_0.01_all'    :  1833,
               'Lr_12.5_cip_F_0.05_all'    :  3828,
               'Lr_12.5_cip_R_0.05_all'    :  3828,
               'Lr_12.5_cip_F_0.01_top'    :  184,
               'Lr_12.5_cip_R_0.01_top'    :  184,
               'Lr_12.5_cip_F_0.01_bot'    :  193,
               'Lr_12.5_cip_R_0.01_bot'    :  193,
                                              
               'Lr_12.5_nocip_F_0.01_all'  :  1833,
               'Lr_12.5_nocip_R_0.01_all'  :  1833,
               'Lr_12.5_nocip_F_0.05_all'  :  3828,
               'Lr_12.5_nocip_R_0.05_all'  :  3828,
               'Lr_12.5_nocip_F_0.01_top'  :  184,
               'Lr_12.5_nocip_R_0.01_top'  :  184,
               'Lr_12.5_nocip_F_0.01_bot'  :  193,
               'Lr_12.5_nocip_R_0.01_bot'  :  193,
                                              
               'Se_5.5_cip_F_0.01_all'     :  1833,
               'Se_5.5_cip_R_0.01_all'     :  1833,
               'Se_5.5_cip_F_0.05_all'     :  3828,
               'Se_5.5_cip_R_0.05_all'     :  3828,
               'Se_5.5_cip_F_0.01_top'     :  184,
               'Se_5.5_cip_R_0.01_top'     :  184,
               'Se_5.5_cip_F_0.01_bot'     :  193,
               'Se_5.5_cip_R_0.01_bot'     :  193,
                                              
               'Se_5.5_nocip_F_0.01_all'   :  1833,
               'Se_5.5_nocip_R_0.01_all'   :  1833,
               'Se_5.5_nocip_F_0.05_all'   :  3828,
               'Se_5.5_nocip_R_0.05_all'   :  3828,
               'Se_5.5_nocip_F_0.01_top'   :  184,
               'Se_5.5_nocip_R_0.01_top'   :  184,
               'Se_5.5_nocip_F_0.01_bot'   :  193,
               'Se_5.5_nocip_R_0.01_bot'   :  193,
                                              
               'Se_12.5_cip_F_0.01_all'    :  1833,
               'Se_12.5_cip_R_0.01_all'    :  1833,
               'Se_12.5_cip_F_0.05_all'    :  3828,
               'Se_12.5_cip_R_0.05_all'    :  3828,
               'Se_12.5_cip_F_0.01_top'    :  184,
               'Se_12.5_cip_R_0.01_top'    :  184,
               'Se_12.5_cip_F_0.01_bot'    :  193,
               'Se_12.5_cip_R_0.01_bot'    :  193,
                                              
               'Se_12.5_nocip_F_0.01_all'  :  1833,
               'Se_12.5_nocip_R_0.01_all'  :  1833,
               'Se_12.5_nocip_F_0.05_all'  :  3828,
               'Se_12.5_nocip_R_0.05_all'  :  3828,
               'Se_12.5_nocip_F_0.01_top'  :  184,
               'Se_12.5_nocip_R_0.01_top'  :  184,
               'Se_12.5_nocip_F_0.01_bot'  :  193,
               'Se_12.5_nocip_R_0.01_bot'  :  193,                                                                                       
               }    
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=4
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(positions[0], positions[-1])
    
    #Plot FE over genes.
    plt.figure(figsize=(4, 4), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    
    ##Transcription units below.
    if set_name=="Lr_5.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Lr_5.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Lr_12.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Se_5.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs['Se_12.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
    elif set_name=="Lr_5.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Lr_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Lr_5.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Lr_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Lr_12.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Se_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Se_5.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions, np.array(dict_of_wigs['Se_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs['Se_12.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
                  
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1, 500).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)   
    plot1.set_xlim([-1500, 1500])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_RPKM_over_{set_type}_FR_{win_width}bp_with_IGR_{length}_bp.png', dpi=400, figsize=(4, 4), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
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
    if set_name=="Lr_5.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_5.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_12.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_5.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5_cip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_cip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_12.5_cip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
    elif set_name=="Lr_5.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.01 ({GSCs_sets["Lr_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 5.5 all bot GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 5.5 all top GCSs ({GSCs_sets["Lr_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_5.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 5.5 all GCSs 0.05 ({GSCs_sets["Lr_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Lr_12.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.01 ({GSCs_sets["Lr_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Lr 12.5 all bot GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Lr 12.5 all top GCSs ({GSCs_sets["Lr_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Lr_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Lr_12.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Lr 12.5 all GCSs 0.05 ({GSCs_sets["Lr_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_5.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.01 ({GSCs_sets["Se_5.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 5.5 all bot GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 5.5 all top GCSs ({GSCs_sets["Se_5.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_5.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_5.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 5.5 all GCSs 0.05 ({GSCs_sets["Se_5.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42      
    elif set_name=="Se_12.5_nocip":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_all']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_all']), linestyle='-', color='#757d8b', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.01 ({GSCs_sets["Se_12.5_cip_F_0.01_all"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_bot']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_bot']), linestyle='-', color='#333738', linewidth=2, alpha=1, label=f'Se 12.5 all bot GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_bot"]})', zorder=2) #R123 -0.1; R12 -0; R3 -0.45             
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.01_top']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.01_top']), linestyle='-', color='#b08642', linewidth=2, alpha=1, label=f'Se 12.5 all top GCSs ({GSCs_sets["Se_12.5_cip_F_0.01_top"]})', zorder=3) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Se_12.5_nocip_F_0.05_all']) + np.array(dict_of_wigs_sm['Se_12.5_nocip_R_0.05_all']), linestyle='-', color='#e4d1b4', linewidth=2, alpha=1, label=f'Se 12.5 all GCSs 0.05 ({GSCs_sets["Se_12.5_cip_F_0.05_all"]})', zorder=4) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42              
      
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1, 500).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.axvline(0, color='black', linestyle=':', alpha=0.7, linewidth=1.5)
    plot1.axvline(length, color='black', linestyle=':', alpha=0.7, linewidth=1.5)    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    plot1.legend(fontsize=11, frameon=False, markerscale=2, handlelength=0.5, handletextpad=0.2, loc='upper left')    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20) 
    plot1.set_xlim([-1500, 1500])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\{set_name}_RPKM_over_{set_type}_FR_smoothed_{win_width}bp_with_IGR_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(4, 4), transparent=True)   
    plt.close()  
    
    return

#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_5.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_5.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_12.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_12.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_5.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_5.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_12.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)
#Name of the signal to plotted (protein or smth.).
Signal_name='Se_12.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict_GSCs, Sm_window, Out_path, Signal_name, Set_type)