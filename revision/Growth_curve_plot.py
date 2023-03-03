###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#Plots growth curve data with standard error intervals.
###############################################

#######
#Packages to be imported.
#######

import random as rd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import scipy
from scipy.stats import norm
import pandas as pd


#########
## Import data for E. coli BL21 harbouring: pBAD30, pBAD30 SeAgo, pBAD30 LrAgo, pBAD30 -Ara.
#########

#Path to the raw data.
Growth_curves_data="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Cfx_titration_growth_curves\Growth_curve_data.xlsx"
#Prefix name of a worksheet.
Prefix_WS_name="Replicate_"
#Path to the output plots.
Outpath_plot="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Cfx_titration_growth_curves\Growth_curve_data"


def prepare_data(replicates):
    mean_ar=[]
    conf_up_ar=[]
    conf_dn_ar=[]
    #Calculate mean and 0.95 confidencial interval.
    for i in range(len(replicates[0])):
        data_point_ar=[]
        for j in range(len(replicates)):
            data_point_ar.append(replicates[j][i])
        data_point_mean=np.mean(data_point_ar)
        data_point_CI=np.std(data_point_ar)*1.96/np.sqrt(len(data_point_ar))
        data_point_up=data_point_mean+data_point_CI
        data_point_dn=data_point_mean-data_point_CI
        mean_ar.append(data_point_mean)
        conf_up_ar.append(data_point_up)
        conf_dn_ar.append(data_point_dn)
           
    return mean_ar, conf_up_ar, conf_dn_ar


def Plot_growth_curve(data_inpath, sheetname_prefix, outpath):
    
    GC_all_data=[]
    #Read growth curves data.
    for i in range(1,4,1):
        gc_data=pd.read_excel(data_inpath, sheet_name=sheetname_prefix+str(i), header=0, index_col=0)
        GC_all_data.append(gc_data)
        
    #Get time points data.
    Time_series=GC_all_data[0].index.tolist()
    print(Time_series)
    
    #Get samples names.
    Samples_names=GC_all_data[0].columns.tolist()
    print(Samples_names)  
    
    Dict_of_mean={}
    Dict_of_conf={}
    for sample_name in Samples_names:
        Replicates=[]
        for j in range(len(GC_all_data)):
            print(f'Sample name: {sample_name}; Replicate number: {j+1}')
            Replicates.append(GC_all_data[j].loc[:, sample_name].values.tolist())
        Mean_ar, Conf_up_ar, Conf_dn_ar=prepare_data(Replicates)
        Dict_of_mean[sample_name]=Mean_ar
        Dict_of_conf[sample_name]=[Conf_up_ar, Conf_dn_ar]
    
    Colors_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    Fill_color_list=['#e07a5f', '#3d405b', '#81b29a', '#f2cc8f']
    color='k'
    
    #Plot data.
    Number_of_plots=11  #15 for full, 11 for cropped.
    Number_of_series=4
    fig, plots=plt.subplots(Number_of_series, Number_of_plots, figsize=(Number_of_plots, Number_of_series), dpi=100)
    rgb = plt.colormaps['rainbow']
    print(rgb(0), rgb(0.5), rgb(1))
    #cmap=plt.get_cmap("jet")
    for i in range(Number_of_series):
        for j in range(Number_of_plots):
            cust_color=rgb(float(j)/Number_of_plots)
            Sample_name=Samples_names[(i*Number_of_plots)+j]
            plots[i, j].plot(Time_series, Dict_of_mean[Sample_name], color=cust_color, linewidth=2, label=Sample_name)
            plots[i, j].plot(Time_series, Dict_of_conf[Sample_name][0], color=cust_color, linewidth=1, alpha=0.6)
            plots[i, j].plot(Time_series, Dict_of_conf[Sample_name][1], color=cust_color, linewidth=1, alpha=0.6)
            plots[i, j].fill_between(Time_series, Dict_of_conf[Sample_name][1], Dict_of_conf[Sample_name][0], color=cust_color, alpha=0.3, linewidth=0.2)
            plots[i, j].spines["top"].set_visible(False)
            plots[i, j].spines["right"].set_visible(False)
            plots[i, j].spines["bottom"].set_linewidth(1.0)
            plots[i, j].spines["left"].set_linewidth(1.0)
            if j==0 and i!=Number_of_series-1:
                Sample_name_ar=Sample_name.split(' ')
                if len(Sample_name_ar)==4:
                    Condition_name=Sample_name_ar[0] + '\n' + Sample_name_ar[1]
                else:
                    Condition_name=Sample_name_ar[0] + '\n' + Sample_name_ar[1] + '\n' + Sample_name_ar[2]
                plots[i, j].set_ylabel(Condition_name + '\nOD$_{600}$', size=8)
                plots[i, j].set_xticks([0, 180, 360, 540, 720]) 
                plots[i, j].set_xticklabels([]) 
                plots[i, j].set_yticks([0,1,2])  
                plots[i, j].set_yticklabels([0,1,2], fontsize=8) 
            elif j!=0 and i==Number_of_series-1:
                plots[i, j].set_xlabel('time, h', size=8)
                plots[i, j].set_xticks([0, 180, 360, 540, 720])
                plots[i, j].set_xticklabels([0, 3, 6, 9, 12], fontsize=8)
                plots[i, j].set_yticks([0,1,2]) 
                plots[i, j].set_yticklabels([])    
            elif j==0 and i==Number_of_series-1:
                Sample_name_ar=Sample_name.split(' ')
                Condition_name=Sample_name_ar[0] + '\n' + Sample_name_ar[1] + '\n' + Sample_name_ar[2]               
                plots[i, j].set_ylabel(Condition_name + '\nOD$_{600}$', size=8)
                plots[i, j].set_xlabel('time, h', size=8)
                plots[i, j].set_xticks([0, 180, 360, 540, 720])
                plots[i, j].set_xticklabels([0, 3, 6, 9, 12], fontsize=8)
                plots[i, j].set_yticks([0,1,2]) 
                plots[i, j].set_yticklabels([0,1,2], fontsize=8)
            else:
                plots[i, j].set_xticks([0, 180, 360, 540, 720]) 
                plots[i, j].set_xticklabels([])  
                plots[i, j].set_yticks([0,1,2]) 
                plots[i, j].set_yticklabels([])    
            plots[i, j].set_ylim([0,2])
            if i==0:
                plots[i, j].set_title('Cfx ' + Sample_name.split(' ')[2] + '\n ng/mL', fontsize=8)
    
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath+'.png', dpi=300, figsize=(Number_of_plots, Number_of_series), transparent=True)
    plt.savefig(outpath+'.svg', dpi=300, figsize=(Number_of_plots, Number_of_series), transparent=True)
    
    return

Plot_growth_curve(Growth_curves_data, Prefix_WS_name, Outpath_plot)
