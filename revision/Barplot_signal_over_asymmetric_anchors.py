###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#Takes files with information about the distribution of some 
#signal over asymmetric anchors (e.g., Chi-sites). Plots this information in form of barplots with statistics.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
import matplotlib.patheffects as PathEffects

#Path to the working directory.
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\'

#Half-window width will be used to smooth signal.
Sm_window=500
#Dictionary of pathes to input WIG data.
Wig_data_in_dict={'Cfx_0.01_Lr_5.5_cip_FN'      :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_cip_RN'      :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_FN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_RN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_FN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_RN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_FN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_RN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_FN'      :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_RN'      :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_FN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_RN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_FN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_RN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_FN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_RN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RN_width_5000bp.txt',
                                                                                                                                                                      
                  'Cfx_0.01_Lr_5.5_cip_FnN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_cip_RnN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_FnN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_RnN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_FnN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_RnN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_FnN'  :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_RnN'  :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_FnN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_RnN'     :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_FnN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_RnN'    :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_FnN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_RnN'   :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_FnN'  :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FnN_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_RnN'  :  PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RnN_width_5000bp.txt',
                                                                                                                                                                     
                  'Cfx_0.01_Lr_5.5_cip_FCUS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_cip_RCUS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_FCUS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_RCUS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_FCUS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_RCUS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_FCUS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_RCUS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_FCUS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_RCUS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_FCUS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_RCUS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_FCUS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_RCUS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_FCUS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FCUS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_RCUS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RCUS_width_5000bp.txt',
                                                                              
                  'Cfx_0.01_Lr_5.5_cip_FCDS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_cip_RCDS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_FCDS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_cip_RCDS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_FCDS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_5.5_nocip_RCDS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_FCDS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Lr_12.5_nocip_RCDS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_FCDS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_cip_RCDS'    :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_cip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_FCDS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_cip_RCDS'   :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_cip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_FCDS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_5.5_nocip_RCDS'  :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_FCDS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FCDS_width_5000bp.txt',
                  'Cfx_0.01_Se_12.5_nocip_RCDS' :  PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig_full\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RCDS_width_5000bp.txt',
                  }

#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="Chi_sites"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Plot_combinations\\'
Dir_check_create(Out_path)


#######
#Parses multyWIG file with RPKM over anchors.
#######

def wig_FE_full_over_genes_parsing(name, full_wigfile):
    print('Now is processing: ' + str(name) + ' ' + str(full_wigfile))
    fullwigin=open(full_wigfile, 'r')
    NE_values=[]
    for line in fullwigin:
        line=line.rstrip(', \n').split('\t')
        Site_name=line[0]
        Site_data=line[1].split(', ')
        for i in range(len(Site_data)):
            Site_data[i]=float(Site_data[i])
        NE_values.append(Site_data)
    fullwigin.close()
    
    return NE_values

#######
#Calculate mean FE, FE STD, and return FE array for further statistics.
#######

def calc_FE_STD(data_cip_F, data_cip_R, data_nocip_F, data_nocip_R, set_name):
    FE_ar_values=[]
    data_cip=data_cip_F+data_cip_R
    data_nocip=data_nocip_F+data_nocip_R
    if len(data_cip_F)+len(data_cip_R)==len(data_nocip_F)+len(data_nocip_R):
        print(f'Number of Chi sites considered: {len(data_cip_F)+len(data_cip_R)}')
    else:
        print(f'Number of Chi sites in cip sample: {len(data_cip_F)+len(data_cip_R)}\nNumber of Chi sites in nocip sample: {len(data_nocip_F)+len(data_nocip_R)}\n')
        return
    for i in range(len(data_cip_F)):
        Site_mean_cip=np.mean(data_cip_F[i][5000+4:])
        Site_mean_nocip=np.mean(data_nocip_F[i][5000+4:])
        if Site_mean_nocip>0:
            Site_FE=Site_mean_cip/Site_mean_nocip
        else:
            Site_FE=(Site_mean_cip+1)/(Site_mean_nocip+1)
        FE_ar_values.append(Site_FE)
        
    for i in range(len(data_cip_R)):
        data_cip_R_i=data_cip_R[i][::-1]
        data_nocip_R_i=data_nocip_R[i][::-1]
        Site_mean_cip=np.mean(data_cip_R_i[5000+4:])
        Site_mean_nocip=np.mean(data_nocip_R_i[5000+4:])
        if Site_mean_nocip>0:
            Site_FE=Site_mean_cip/Site_mean_nocip
        else:
            Site_FE=(Site_mean_cip+1)/(Site_mean_nocip+1)
        FE_ar_values.append(Site_FE)
    
    FE_mean_value=np.mean(FE_ar_values)
    FE_STD_value=np.std(FE_ar_values)
    
    print(f'For sample {set_name}, FE mean={round(FE_mean_value, 2)}, FE STD={round(FE_STD_value, 2)}')
    
    return FE_mean_value, FE_STD_value, FE_ar_values


#######
#Plot the signal for different groups of anchors together.
#Fold enrichment or other ratios.
#######


def plot_FE_over_anchors(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    Chi_sets={'Cfx_0.01_Lr_5.5_cip_FN'      :  198,
              'Cfx_0.01_Lr_5.5_cip_RN'      :  188,
              'Cfx_0.01_Lr_12.5_cip_FN'     :  198,
              'Cfx_0.01_Lr_12.5_cip_RN'     :  188,
              'Cfx_0.01_Lr_5.5_nocip_FN'    :  198,
              'Cfx_0.01_Lr_5.5_nocip_RN'    :  188,
              'Cfx_0.01_Lr_12.5_nocip_FN'   :  198,
              'Cfx_0.01_Lr_12.5_nocip_RN'   :  188,
              'Cfx_0.01_Se_5.5_cip_FN'      :  198,
              'Cfx_0.01_Se_5.5_cip_RN'      :  188,
              'Cfx_0.01_Se_12.5_cip_FN'     :  198,
              'Cfx_0.01_Se_12.5_cip_RN'     :  188,
              'Cfx_0.01_Se_5.5_nocip_FN'    :  198,
              'Cfx_0.01_Se_5.5_nocip_RN'    :  188,
              'Cfx_0.01_Se_12.5_nocip_FN'   :  198,
              'Cfx_0.01_Se_12.5_nocip_RN'   :  188,
                                            
              'Cfx_0.01_Lr_5.5_cip_FnN'     :  197,
              'Cfx_0.01_Lr_5.5_cip_RnN'     :  250,
              'Cfx_0.01_Lr_12.5_cip_FnN'    :  197,
              'Cfx_0.01_Lr_12.5_cip_RnN'    :  250,
              'Cfx_0.01_Lr_5.5_nocip_FnN'   :  197,
              'Cfx_0.01_Lr_5.5_nocip_RnN'   :  250,
              'Cfx_0.01_Lr_12.5_nocip_FnN'  :  197,
              'Cfx_0.01_Lr_12.5_nocip_RnN'  :  250,
              'Cfx_0.01_Se_5.5_cip_FnN'     :  197,
              'Cfx_0.01_Se_5.5_cip_RnN'     :  250,
              'Cfx_0.01_Se_12.5_cip_FnN'    :  197,
              'Cfx_0.01_Se_12.5_cip_RnN'    :  250,
              'Cfx_0.01_Se_5.5_nocip_FnN'   :  197,
              'Cfx_0.01_Se_5.5_nocip_RnN'   :  250,
              'Cfx_0.01_Se_12.5_nocip_FnN'  :  197,
              'Cfx_0.01_Se_12.5_nocip_RnN'  :  250,
              
              'Cfx_0.01_Lr_5.5_cip_FCUS'    :  90,
              'Cfx_0.01_Lr_5.5_cip_RCUS'    :  101,
              'Cfx_0.01_Lr_12.5_cip_FCUS'   :  90,
              'Cfx_0.01_Lr_12.5_cip_RCUS'   :  101,
              'Cfx_0.01_Lr_5.5_nocip_FCUS'  :  90,
              'Cfx_0.01_Lr_5.5_nocip_RCUS'  :  101,
              'Cfx_0.01_Lr_12.5_nocip_FCUS' :  90,
              'Cfx_0.01_Lr_12.5_nocip_RCUS' :  101,
              'Cfx_0.01_Se_5.5_cip_FCUS'    :  90,
              'Cfx_0.01_Se_5.5_cip_RCUS'    :  101,
              'Cfx_0.01_Se_12.5_cip_FCUS'   :  90,
              'Cfx_0.01_Se_12.5_cip_RCUS'   :  101,
              'Cfx_0.01_Se_5.5_nocip_FCUS'  :  90,
              'Cfx_0.01_Se_5.5_nocip_RCUS'  :  101,
              'Cfx_0.01_Se_12.5_nocip_FCUS' :  90,
              'Cfx_0.01_Se_12.5_nocip_RCUS' :  101,
                                               
              'Cfx_0.01_Lr_5.5_cip_FCDS'    :  90,
              'Cfx_0.01_Lr_5.5_cip_RCDS'    :  102,
              'Cfx_0.01_Lr_12.5_cip_FCDS'   :  90,
              'Cfx_0.01_Lr_12.5_cip_RCDS'   :  102,
              'Cfx_0.01_Lr_5.5_nocip_FCDS'  :  90,
              'Cfx_0.01_Lr_5.5_nocip_RCDS'  :  102,
              'Cfx_0.01_Lr_12.5_nocip_FCDS' :  90,
              'Cfx_0.01_Lr_12.5_nocip_RCDS' :  102,
              'Cfx_0.01_Se_5.5_cip_FCDS'    :  90,
              'Cfx_0.01_Se_5.5_cip_RCDS'    :  102,
              'Cfx_0.01_Se_12.5_cip_FCDS'   :  90,
              'Cfx_0.01_Se_12.5_cip_RCDS'   :  102,
              'Cfx_0.01_Se_5.5_nocip_FCDS'  :  90,
              'Cfx_0.01_Se_5.5_nocip_RCDS'  :  102,
              'Cfx_0.01_Se_12.5_nocip_FCDS' :  90,
              'Cfx_0.01_Se_12.5_nocip_RCDS' :  102,                    
              }
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    for name, file in wig_in_dict.items():
        data=wig_FE_full_over_genes_parsing(name, file)
        dict_of_wigs[name]=data                  
    
    ##Anchor sets below.
    Dataset_name='Lr_5.5'
    Lr_5_5_FE_mean_value_FN,   Lr_5_5_FE_STD_value_FN,   Lr_5_5_FE_ar_values_FN  =calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FN'],   dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RN'],   dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FN'],   dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RN'],   'Lr_5.5_N')
    Lr_5_5_FE_mean_value_FnN,  Lr_5_5_FE_STD_value_FnN,  Lr_5_5_FE_ar_values_FnN =calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FnN'],  dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RnN'],  dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FnN'],  dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RnN'],  'Lr_5.5_nN')
    Lr_5_5_FE_mean_value_FCUS, Lr_5_5_FE_STD_value_FCUS, Lr_5_5_FE_ar_values_FCUS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FCUS'], dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RCUS'], dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FCUS'], dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RCUS'], 'Lr_5.5_US')
    Lr_5_5_FE_mean_value_FCDS, Lr_5_5_FE_STD_value_FCDS, Lr_5_5_FE_ar_values_FCDS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FCDS'], dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RCDS'], dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FCDS'], dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RCDS'], 'Lr_5.5_DS')
    
    Intervals_stat=stats.ttest_ind(Lr_5_5_FE_ar_values_FN, Lr_5_5_FE_ar_values_FnN, equal_var=True)
    print(f'Working with dataset: {Dataset_name}')
    print('Comparing of Chi sites enrichment with adjacent GCS vs Chi sites enrichment with no adjacent GCSs')
    print(f'Sample size: {len(Lr_5_5_FE_ar_values_FN)}, Sample size: {len(Lr_5_5_FE_ar_values_FnN)}')
    print(f'\nT-test FE, Mean1={round(Lr_5_5_FE_mean_value_FN,3)}; Mean2={round(Lr_5_5_FE_mean_value_FnN,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    Intervals_stat=stats.ttest_ind(Lr_5_5_FE_ar_values_FCUS, Lr_5_5_FE_ar_values_FCDS, equal_var=True)
    print('Comparing of Chi sites enrichment with a GCS upstream vs Chi sites enrichment with GCS downstream')
    print(f'Sample size: {len(Lr_5_5_FE_ar_values_FCUS)}, Sample size: {len(Lr_5_5_FE_ar_values_FCDS)}')
    print(f'\nT-test FE, Mean1={round(Lr_5_5_FE_mean_value_FCUS,3)}; Mean2={round(Lr_5_5_FE_mean_value_FCDS,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    Dataset_name='Lr_12.5'
    Lr_12_5_FE_mean_value_FN,   Lr_12_5_FE_STD_value_FN,   Lr_12_5_FE_ar_values_FN  =calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FN'],   dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RN'],   dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FN'],   dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RN'],   'Lr_12.5_N')
    Lr_12_5_FE_mean_value_FnN,  Lr_12_5_FE_STD_value_FnN,  Lr_12_5_FE_ar_values_FnN =calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FnN'],  dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RnN'],  dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FnN'],  dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RnN'],  'Lr_12.5_nN')
    Lr_12_5_FE_mean_value_FCUS, Lr_12_5_FE_STD_value_FCUS, Lr_12_5_FE_ar_values_FCUS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FCUS'], dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RCUS'], dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FCUS'], dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RCUS'], 'Lr_12.5_US')
    Lr_12_5_FE_mean_value_FCDS, Lr_12_5_FE_STD_value_FCDS, Lr_12_5_FE_ar_values_FCDS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FCDS'], dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RCDS'], dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FCDS'], dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RCDS'], 'Lr_12.5_DS')
    
    Intervals_stat=stats.ttest_ind(Lr_12_5_FE_ar_values_FN, Lr_12_5_FE_ar_values_FnN, equal_var=True)
    print(f'Working with dataset: {Dataset_name}')
    print('Comparing of Chi sites enrichment with adjacent GCS vs Chi sites enrichment with no adjacent GCSs')
    print(f'Sample size: {len(Lr_12_5_FE_ar_values_FN)}, Sample size: {len(Lr_12_5_FE_ar_values_FnN)}')
    print(f'\nT-test FE, Mean1={round(Lr_12_5_FE_mean_value_FN,3)}; Mean2={round(Lr_12_5_FE_mean_value_FnN,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    Intervals_stat=stats.ttest_ind(Lr_12_5_FE_ar_values_FCUS, Lr_12_5_FE_ar_values_FCDS, equal_var=True)
    print('Comparing of Chi sites enrichment with a GCS upstream vs Chi sites enrichment with GCS downstream')
    print(f'Sample size: {len(Lr_12_5_FE_ar_values_FCUS)}, Sample size: {len(Lr_12_5_FE_ar_values_FCDS)}')
    print(f'\nT-test FE, Mean1={round(Lr_12_5_FE_mean_value_FCUS,3)}; Mean2={round(Lr_12_5_FE_mean_value_FCDS,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
     
    Dataset_name='Se_5.5'
    Se_5_5_FE_mean_value_FN,   Se_5_5_FE_STD_value_FN,   Se_5_5_FE_ar_values_FN  =calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FN'],   dict_of_wigs['Cfx_0.01_Se_5.5_cip_RN'],   dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FN'],   dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RN'],   'Se_5.5_N')
    Se_5_5_FE_mean_value_FnN,  Se_5_5_FE_STD_value_FnN,  Se_5_5_FE_ar_values_FnN =calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FnN'],  dict_of_wigs['Cfx_0.01_Se_5.5_cip_RnN'],  dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FnN'],  dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RnN'],  'Se_5.5_nN')
    Se_5_5_FE_mean_value_FCUS, Se_5_5_FE_STD_value_FCUS, Se_5_5_FE_ar_values_FCUS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FCUS'], dict_of_wigs['Cfx_0.01_Se_5.5_cip_RCUS'], dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FCUS'], dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RCUS'], 'Se_5.5_US')
    Se_5_5_FE_mean_value_FCDS, Se_5_5_FE_STD_value_FCDS, Se_5_5_FE_ar_values_FCDS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FCDS'], dict_of_wigs['Cfx_0.01_Se_5.5_cip_RCDS'], dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FCDS'], dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RCDS'], 'Se_5.5_DS')
    
    Intervals_stat=stats.ttest_ind(Se_5_5_FE_ar_values_FN, Se_5_5_FE_ar_values_FnN, equal_var=True)
    print(f'Working with dataset: {Dataset_name}')
    print('Comparing of Chi sites enrichment with adjacent GCS vs Chi sites enrichment with no adjacent GCSs')
    print(f'Sample size: {len(Se_5_5_FE_ar_values_FN)}, Sample size: {len(Se_5_5_FE_ar_values_FnN)}')
    print(f'\nT-test FE, Mean1={round(Se_5_5_FE_mean_value_FN,3)}; Mean2={round(Se_5_5_FE_mean_value_FnN,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    Intervals_stat=stats.ttest_ind(Se_5_5_FE_ar_values_FCUS, Se_5_5_FE_ar_values_FCDS, equal_var=True)
    print('Comparing of Chi sites enrichment with a GCS upstream vs Chi sites enrichment with GCS downstream')
    print(f'Sample size: {len(Se_5_5_FE_ar_values_FCUS)}, Sample size: {len(Se_5_5_FE_ar_values_FCDS)}')
    print(f'\nT-test FE, Mean1={round(Se_5_5_FE_mean_value_FCUS,3)}; Mean2={round(Se_5_5_FE_mean_value_FCDS,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
       
    Dataset_name='Se_12.5'
    Se_12_5_FE_mean_value_FN,   Se_12_5_FE_STD_value_FN,   Se_12_5_FE_ar_values_FN  =calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FN'],   dict_of_wigs['Cfx_0.01_Se_12.5_cip_RN'],   dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FN'],   dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RN'],   'Se_12.5_N')
    Se_12_5_FE_mean_value_FnN,  Se_12_5_FE_STD_value_FnN,  Se_12_5_FE_ar_values_FnN =calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FnN'],  dict_of_wigs['Cfx_0.01_Se_12.5_cip_RnN'],  dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FnN'],  dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RnN'],  'Se_12.5_nN')
    Se_12_5_FE_mean_value_FCUS, Se_12_5_FE_STD_value_FCUS, Se_12_5_FE_ar_values_FCUS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FCUS'], dict_of_wigs['Cfx_0.01_Se_12.5_cip_RCUS'], dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FCUS'], dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RCUS'], 'Se_12.5_US')
    Se_12_5_FE_mean_value_FCDS, Se_12_5_FE_STD_value_FCDS, Se_12_5_FE_ar_values_FCDS=calc_FE_STD(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FCDS'], dict_of_wigs['Cfx_0.01_Se_12.5_cip_RCDS'], dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FCDS'], dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RCDS'], 'Se_12.5_DS')
    
    Intervals_stat=stats.ttest_ind(Se_12_5_FE_ar_values_FN, Se_12_5_FE_ar_values_FnN, equal_var=True)
    print(f'Working with dataset: {Dataset_name}')
    print('Comparing of Chi sites enrichment with adjacent GCS vs Chi sites enrichment with no adjacent GCSs')
    print(f'Sample size: {len(Se_12_5_FE_ar_values_FN)}, Sample size: {len(Se_12_5_FE_ar_values_FnN)}')
    print(f'\nT-test FE, Mean1={round(Se_12_5_FE_mean_value_FN,3)}; Mean2={round(Se_12_5_FE_mean_value_FnN,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
    
    Intervals_stat=stats.ttest_ind(Se_12_5_FE_ar_values_FCUS, Se_12_5_FE_ar_values_FCDS, equal_var=True)
    print('Comparing of Chi sites enrichment with a GCS upstream vs Chi sites enrichment with GCS downstream')
    print(f'Sample size: {len(Se_12_5_FE_ar_values_FCUS)}, Sample size: {len(Se_12_5_FE_ar_values_FCDS)}')
    print(f'\nT-test FE, Mean1={round(Se_12_5_FE_mean_value_FCUS,3)}; Mean2={round(Se_12_5_FE_mean_value_FCDS,3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 
                           
    #Compare anchor sets.
    #Plot FE over anchors. GCS near vs no GCSs near a Chi site.
    fig, plot_av=plt.subplots(1,1,figsize=(3,3), dpi=100)

    Mean_data_ar=[Lr_5_5_FE_mean_value_FN, Lr_5_5_FE_mean_value_FnN, Lr_12_5_FE_mean_value_FN, Lr_12_5_FE_mean_value_FnN, Se_5_5_FE_mean_value_FN, Se_5_5_FE_mean_value_FnN, Se_12_5_FE_mean_value_FN, Se_12_5_FE_mean_value_FnN]
    STD_data=[Lr_5_5_FE_STD_value_FN, Lr_5_5_FE_STD_value_FnN, Lr_12_5_FE_STD_value_FN, Lr_12_5_FE_STD_value_FnN, Se_5_5_FE_STD_value_FN, Se_5_5_FE_STD_value_FnN, Se_12_5_FE_STD_value_FN, Se_12_5_FE_STD_value_FnN]
    Loc_max_ar=[]
    for i in range(len(Mean_data_ar)):
        Loc_max=Mean_data_ar[i]+STD_data[i]
        Loc_max_ar.append(Loc_max)
    Glob_max=np.max(Loc_max_ar)
    
    Conditions=['near\nGCS', 'no\nGSC', 'near\nGCS', 'no\nGSC', 'near\nGCS', 'no\nGSC', 'near\nGCS', 'no\nGSC']
    
    color_list=['#c96458', '#6a65c7', '#c96458', '#6a65c7', '#c96458', '#6a65c7', '#c96458', '#6a65c7']
    X_coords=[1, 2, 4, 5, 7, 8, 10, 11]
    X_coords_m=[1, 2, 4, 5, 7, 8, 10, 11]
      
    Bars=plot_av.bar(X_coords, Mean_data_ar, yerr=STD_data, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.9, color=color_list, edgecolor='k', linewidth=1)
    plot_av.set_ylabel('RPKM cip/RPKM nocip', size=16)
    plot_av.set_xticks(X_coords_m)
    plot_av.set_xticklabels(Conditions, rotation=0, size=9)     
    plot_av.set_ylim([0.45, Glob_max*1.1])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\Barplot_N_vs_nN_{set_name}_over_{set_type}_cip_div_nocip_FR_5000_bp_US.png', dpi=400, figsize=(3, 3), transparent=True)   
    plt.savefig(f'{output_path}\\Barplot_N_vs_nN_{set_name}_over_{set_type}_cip_div_nocip_FR_5000_bp_US.svg', dpi=400, figsize=(3, 3), transparent=True)   
    plt.close() 
    
    
    #Plot FE over anchors. Chi upstream of a GCS vs Chi downstream of a GCS site.
    fig, plot_av=plt.subplots(1,1,figsize=(3,3), dpi=100)

    Mean_data_ar=[Lr_5_5_FE_mean_value_FCUS, Lr_5_5_FE_mean_value_FCDS, Lr_12_5_FE_mean_value_FCUS, Lr_12_5_FE_mean_value_FCDS, Se_5_5_FE_mean_value_FCUS, Se_5_5_FE_mean_value_FCDS, Se_12_5_FE_mean_value_FCUS, Se_12_5_FE_mean_value_FCDS]
    STD_data=[Lr_5_5_FE_STD_value_FCUS, Lr_5_5_FE_STD_value_FCDS, Lr_12_5_FE_STD_value_FCUS, Lr_12_5_FE_STD_value_FCDS, Se_5_5_FE_STD_value_FCUS, Se_5_5_FE_STD_value_FCDS, Se_12_5_FE_STD_value_FCUS, Se_12_5_FE_STD_value_FCDS]
    Loc_max_ar=[]
    for i in range(len(Mean_data_ar)):
        Loc_max=Mean_data_ar[i]+STD_data[i]
        Loc_max_ar.append(Loc_max)
    Glob_max=np.max(Loc_max_ar)    
    
    Conditions=['GCS\nUS', 'GCS\nDS', 'GCS\nUS', 'GCS\nDS', 'GCS\nUS', 'GCS\nDS', 'GCS\nUS', 'GCS\nDS']
    
    color_list=['#c96458', '#6a65c7', '#c96458', '#6a65c7', '#c96458', '#6a65c7', '#c96458', '#6a65c7']
    X_coords=[1, 2, 4, 5, 7, 8, 10, 11]
    X_coords_m=[1, 2, 4, 5, 7, 8, 10, 11]
      
    Bars=plot_av.bar(X_coords, Mean_data_ar, yerr=STD_data, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.9, color=color_list, edgecolor='k', linewidth=1)
    plot_av.set_ylabel('RPKM cip/RPKM nocip', size=16)
    plot_av.set_xticks(X_coords_m)
    plot_av.set_xticklabels(Conditions, rotation=0, size=9)     
    plot_av.set_ylim([0.45, Glob_max*1.1])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\Barplot_US_vs_DS_{set_name}_over_{set_type}_cip_div_nocip_FR_5000_bp_US.png', dpi=400, figsize=(3, 3), transparent=True)   
    plt.savefig(f'{output_path}\\Barplot_US_vs_DS_{set_name}_over_{set_type}_cip_div_nocip_FR_5000_bp_US.svg', dpi=400, figsize=(3, 3), transparent=True)   
    plt.close() 
    
    return


#Lr Se Cfx effects.
#Name of the signal to plotted (protein or smth.).
Signal_name='All_Ago_Seq_sets'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
