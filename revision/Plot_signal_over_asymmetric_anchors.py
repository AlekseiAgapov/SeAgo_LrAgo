###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#signal over asymmetric anchors (e.g., Chi-sites). Plots this information.

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
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\'

#Half-window width will be used to smooth signal.
Sm_window=1000
#Dictionary of pathes to input WIG data.
Wig_data_in_dict={'Cfx_0.01_Lr_5.5_cip_F'               :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_R'               :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_F'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_R'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_F'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_R'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_F'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_R'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_F'               :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_R'               :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_F'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_R'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_F'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_R'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_F'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_R'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_All_R_width_50000bp.wig',
                  
                  'Cfx_0.01_Lr_5.5_cip_FR'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_RF'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_FR'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_RF'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_FR'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_RF'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_FR'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_RF'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_FR'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_RF'              :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_FR'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_RF'             :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_FR'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_RF'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_All_F_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_FR'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_All_R_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_RF'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_All_F_width_50000bp.wig',
                  
                  'Cfx_0.01_Lr_5.5_cip_FN'              :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_RN'              :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_FN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_RN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_FN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_RN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_FN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_RN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_FN'              :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_RN'              :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_FN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_RN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_FN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_RN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_FN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_RN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RN_width_50000bp.wig',
                                                                                                                                                                                
                  'Cfx_0.01_Lr_5.5_cip_FnN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_RnN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_FnN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_RnN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_FnN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_RnN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_FnN'          :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_RnN'          :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_FnN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_RnN'             :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_FnN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_RnN'            :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_FnN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_RnN'           :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_FnN'          :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FnN_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_RnN'          :    PWD + 'Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RnN_width_50000bp.wig',
                                                                                                                                                                                
                  'Cfx_0.01_Lr_5.5_cip_FCUS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_RCUS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_FCUS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_RCUS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_FCUS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_RCUS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_FCUS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_RCUS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_FCUS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_RCUS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_FCUS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_RCUS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_FCUS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_RCUS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_FCUS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FCUS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_RCUS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RCUS_width_50000bp.wig',
                                                                                        
                  'Cfx_0.01_Lr_5.5_cip_FCDS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_RCDS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_FCDS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_RCDS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_FCDS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_RCDS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_FCDS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_RCDS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_FCDS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_RCDS'            :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_FCDS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_RCDS'           :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_FCDS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_RCDS'          :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_FCDS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_F_over_Cfx_0.01_FCDS_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_RCDS'         :    PWD + 'Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_R_over_Cfx_0.01_RCDS_width_50000bp.wig',

                  'Cfx_0.01_Lr_5.5_cip_US'            :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_cip_OTH'           :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_cip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_US'           :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_cip_OTH'          :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_cip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_US'          :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Lr_5.5_nocip_OTH'         :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_5.5_nocip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_US'         :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Lr_12.5_nocip_OTH'        :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Lr_12.5_nocip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_US'            :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_cip_OTH'           :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_cip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_US'           :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_cip_OTH'          :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_cip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_US'          :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Se_5.5_nocip_OTH'         :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_5.5_nocip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_US'         :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_over_Cfx_0.01_US_width_50000bp.wig',
                  'Cfx_0.01_Se_12.5_nocip_OTH'        :    PWD + 'Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Signal_of_GCSs_wig\Signal_Se_12.5_nocip_over_Cfx_0.01_OTH_width_50000bp.wig',
                  

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
            ww_l=line[2].split('=')[1].rstrip('"').lstrip('"')
            win_width=int(ww_l)
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
#Plot the signal for different groups of anchors together.
#Fold enrichment or other ratios.
#######


def plot_FE_over_anchors(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    Chi_sets={'Cfx_0.01_Lr_5.5_cip_F'    :    395,
              'Cfx_0.01_Lr_5.5_cip_R'    :    438,
              'Cfx_0.01_Lr_12.5_cip_F'   :    395,
              'Cfx_0.01_Lr_12.5_cip_R'   :    438,
              'Cfx_0.01_Lr_5.5_nocip_F'  :    395,
              'Cfx_0.01_Lr_5.5_nocip_R'  :    438,
              'Cfx_0.01_Lr_12.5_nocip_F' :    395,
              'Cfx_0.01_Lr_12.5_nocip_R' :    438,
              'Cfx_0.01_Se_5.5_cip_F'    :    395,
              'Cfx_0.01_Se_5.5_cip_R'    :    438,
              'Cfx_0.01_Se_12.5_cip_F'   :    395,
              'Cfx_0.01_Se_12.5_cip_R'   :    438,
              'Cfx_0.01_Se_5.5_nocip_F'  :    395,
              'Cfx_0.01_Se_5.5_nocip_R'  :    438,
              'Cfx_0.01_Se_12.5_nocip_F' :    395,
              'Cfx_0.01_Se_12.5_nocip_R' :    438,
              'Cfx_0.01_Lr_5.5_cip_FR'    :   395,
              'Cfx_0.01_Lr_5.5_cip_RF'    :   438,
              'Cfx_0.01_Lr_12.5_cip_FR'   :   395,
              'Cfx_0.01_Lr_12.5_cip_RF'   :   438,
              'Cfx_0.01_Lr_5.5_nocip_FR'  :   395,
              'Cfx_0.01_Lr_5.5_nocip_RF'  :   438,
              'Cfx_0.01_Lr_12.5_nocip_FR' :   395,
              'Cfx_0.01_Lr_12.5_nocip_RF' :   438,
              'Cfx_0.01_Se_5.5_cip_FR'    :   395,
              'Cfx_0.01_Se_5.5_cip_RF'    :   438,
              'Cfx_0.01_Se_12.5_cip_FR'   :   395,
              'Cfx_0.01_Se_12.5_cip_RF'   :   438,
              'Cfx_0.01_Se_5.5_nocip_FR'  :   395,
              'Cfx_0.01_Se_5.5_nocip_RF'  :   438,
              'Cfx_0.01_Se_12.5_nocip_FR' :   395,
              'Cfx_0.01_Se_12.5_nocip_RF' :   438,
              'Cfx_0.01_Lr_5.5_cip_FN'    :    198,
              'Cfx_0.01_Lr_5.5_cip_RN'    :    188,
              'Cfx_0.01_Lr_12.5_cip_FN'   :    198,
              'Cfx_0.01_Lr_12.5_cip_RN'   :    188,
              'Cfx_0.01_Lr_5.5_nocip_FN'  :    198,
              'Cfx_0.01_Lr_5.5_nocip_RN'  :    188,
              'Cfx_0.01_Lr_12.5_nocip_FN' :    198,
              'Cfx_0.01_Lr_12.5_nocip_RN' :    188,
              'Cfx_0.01_Se_5.5_cip_FN'    :    198,
              'Cfx_0.01_Se_5.5_cip_RN'    :    188,
              'Cfx_0.01_Se_12.5_cip_FN'   :    198,
              'Cfx_0.01_Se_12.5_cip_RN'   :    188,
              'Cfx_0.01_Se_5.5_nocip_FN'  :    198,
              'Cfx_0.01_Se_5.5_nocip_RN'  :    188,
              'Cfx_0.01_Se_12.5_nocip_FN' :    198,
              'Cfx_0.01_Se_12.5_nocip_RN' :    188,
                            
              'Cfx_0.01_Lr_5.5_cip_FnN'    :    197,
              'Cfx_0.01_Lr_5.5_cip_RnN'    :    250,
              'Cfx_0.01_Lr_12.5_cip_FnN'   :    197,
              'Cfx_0.01_Lr_12.5_cip_RnN'   :    250,
              'Cfx_0.01_Lr_5.5_nocip_FnN'  :    197,
              'Cfx_0.01_Lr_5.5_nocip_RnN'  :    250,
              'Cfx_0.01_Lr_12.5_nocip_FnN' :    197,
              'Cfx_0.01_Lr_12.5_nocip_RnN' :    250,
              'Cfx_0.01_Se_5.5_cip_FnN'    :    197,
              'Cfx_0.01_Se_5.5_cip_RnN'    :    250,
              'Cfx_0.01_Se_12.5_cip_FnN'   :    197,
              'Cfx_0.01_Se_12.5_cip_RnN'   :    250,
              'Cfx_0.01_Se_5.5_nocip_FnN'  :    197,
              'Cfx_0.01_Se_5.5_nocip_RnN'  :    250,
              'Cfx_0.01_Se_12.5_nocip_FnN' :    197,
              'Cfx_0.01_Se_12.5_nocip_RnN' :    250,
                            
              'Cfx_0.01_Lr_5.5_cip_FCUS'            :    90,
              'Cfx_0.01_Lr_5.5_cip_RCUS'            :    101,
              'Cfx_0.01_Lr_12.5_cip_FCUS'           :    90,
              'Cfx_0.01_Lr_12.5_cip_RCUS'           :    101,
              'Cfx_0.01_Lr_5.5_nocip_FCUS'          :    90,
              'Cfx_0.01_Lr_5.5_nocip_RCUS'          :    101,
              'Cfx_0.01_Lr_12.5_nocip_FCUS'         :    90,
              'Cfx_0.01_Lr_12.5_nocip_RCUS'         :    101,
              'Cfx_0.01_Se_5.5_cip_FCUS'            :    90,
              'Cfx_0.01_Se_5.5_cip_RCUS'            :    101,
              'Cfx_0.01_Se_12.5_cip_FCUS'           :    90,
              'Cfx_0.01_Se_12.5_cip_RCUS'           :    101,
              'Cfx_0.01_Se_5.5_nocip_FCUS'          :    90,
              'Cfx_0.01_Se_5.5_nocip_RCUS'          :    101,
              'Cfx_0.01_Se_12.5_nocip_FCUS'         :    90,
              'Cfx_0.01_Se_12.5_nocip_RCUS'         :    101,
                                                         
              'Cfx_0.01_Lr_5.5_cip_FCDS'            :    90,
              'Cfx_0.01_Lr_5.5_cip_RCDS'            :    102,
              'Cfx_0.01_Lr_12.5_cip_FCDS'           :    90,
              'Cfx_0.01_Lr_12.5_cip_RCDS'           :    102,
              'Cfx_0.01_Lr_5.5_nocip_FCDS'          :    90,
              'Cfx_0.01_Lr_5.5_nocip_RCDS'          :    102,
              'Cfx_0.01_Lr_12.5_nocip_FCDS'         :    90,
              'Cfx_0.01_Lr_12.5_nocip_RCDS'         :    102,
              'Cfx_0.01_Se_5.5_cip_FCDS'            :    90,
              'Cfx_0.01_Se_5.5_cip_RCDS'            :    102,
              'Cfx_0.01_Se_12.5_cip_FCDS'           :    90,
              'Cfx_0.01_Se_12.5_cip_RCDS'           :    102,
              'Cfx_0.01_Se_5.5_nocip_FCDS'          :    90,
              'Cfx_0.01_Se_5.5_nocip_RCDS'          :    102,
              'Cfx_0.01_Se_12.5_nocip_FCDS'         :    90,
              'Cfx_0.01_Se_12.5_nocip_RCDS'         :    102,
              
              'Cfx_0.01_Lr_5.5_cip_US'            :    188,
              'Cfx_0.01_Lr_5.5_cip_OTH'           :    645,
              'Cfx_0.01_Lr_12.5_cip_US'           :    188,
              'Cfx_0.01_Lr_12.5_cip_OTH'          :    645,
              'Cfx_0.01_Lr_5.5_nocip_US'          :    188,
              'Cfx_0.01_Lr_5.5_nocip_OTH'         :    645,
              'Cfx_0.01_Lr_12.5_nocip_US'         :    188,
              'Cfx_0.01_Lr_12.5_nocip_OTH'        :    645,
              'Cfx_0.01_Se_5.5_cip_US'            :    188,
              'Cfx_0.01_Se_5.5_cip_OTH'           :    645,
              'Cfx_0.01_Se_12.5_cip_US'           :    188,
              'Cfx_0.01_Se_12.5_cip_OTH'          :    645,
              'Cfx_0.01_Se_5.5_nocip_US'          :    188,
              'Cfx_0.01_Se_5.5_nocip_OTH'         :    645,
              'Cfx_0.01_Se_12.5_nocip_US'         :    188,
              'Cfx_0.01_Se_12.5_nocip_OTH'        :    645,
                                                 
              }
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=30000
    length=5000
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+4, 1)    
    print(positions[0], positions[-1])
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    ##Anchor sets below.
    if set_name=='Lr_5.5_FR_FE': 
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_R'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Lr_12.5_FR_FE':    
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_R'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Se_5.5_FR_FE':    
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_R'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Se_12.5_FR_FE':    
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_R'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07            
    elif set_name=='Lr_5.5_FR_vicinity_FE':     
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RN'][::-1])) / (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RN'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR near GCS ({Chi_sets["Cfx_0.01_Lr_5.5_cip_FN"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RnN'][::-1])),    linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR not near GCS ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_12.5_FR_vicinity_FE':   
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RN'][::-1])) / (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR near Chi ({Chi_sets["Cfx_0.01_Lr_12.5_cip_FN"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR not near Chi ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_12.5_FR_vicinity_FE':    
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_RN'][::-1])) / (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR near Chi ({Chi_sets["Cfx_0.01_Se_12.5_cip_FN"]+Chi_sets["Cfx_0.01_Se_12.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR not near Chi ({Chi_sets["Cfx_0.01_Se_12.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_5.5_FR_vicinity_FE':    
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_RN'][::-1])) / (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FN']) + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR near GCS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FN"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FnN'])+ np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR not near GCS ({Chi_sets["Cfx_0.01_Se_5.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_5.5_FR_Cl_FE':      
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_R'][::-1]))   /(np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Lr_12.5_FR_Cl_FE':       
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_R'][::-1]))   /(np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_5.5_FR_Cl_FE':      
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_F'])  + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_R'][::-1]))   /(np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_F'])  + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_12.5_FR_Cl_FE':   
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FCUS'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FCDS'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions, (np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_R'][::-1]))   /(np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_F'])   +np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Lr_5.5_US_vs_OTH':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_US']) /np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip US ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_OTH'])/np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_US']) /np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip US ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_OTH'])/np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_US']) /np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip US ({Chi_sets["Cfx_0.01_Se_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_OTH'])/np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Se_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_US']) /np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip US ({Chi_sets["Cfx_0.01_Se_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_OTH'])/np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Se_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    
    elif set_name=='Lr_5.5_US_vs_OTH_nocip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip DS ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH_nocip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip DS ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH_nocip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip DS ({Chi_sets["Cfx_0.01_Se_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip OTH ({Chi_sets["Cfx_0.01_Se_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH_nocip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip DS ({Chi_sets["Cfx_0.01_Se_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip OTH ({Chi_sets["Cfx_0.01_Se_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    
    elif set_name=='Lr_5.5_US_vs_OTH_cip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip DS ({Chi_sets["Cfx_0.01_Lr_5.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH_cip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip DS ({Chi_sets["Cfx_0.01_Lr_12.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH_cip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip OTH ({Chi_sets["Cfx_0.01_Se_5.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH_cip':      
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip DS ({Chi_sets["Cfx_0.01_Se_12.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip OTH ({Chi_sets["Cfx_0.01_Se_12.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
           
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+4+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    plot1.tick_params(axis='both', labelsize=20)
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_xlim([-30000, 30008])
    if set_name in ['Se_12.5_US_vs_OTH_cip', 'Se_5.5_US_vs_OTH_cip', 'Lr_12.5_US_vs_OTH_cip', 'Lr_5.5_US_vs_OTH_cip']:
        plot1.set_ylabel(f'RPKM cip', size=20)
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_FR_{win_width}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_FR_{win_width}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6                
    elif set_name in ['Se_12.5_US_vs_OTH_nocip', 'Se_5.5_US_vs_OTH_nocip', 'Lr_12.5_US_vs_OTH_nocip', 'Lr_5.5_US_vs_OTH_nocip']:
        plot1.set_ylabel(f'RPKM nocip', size=20)  
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_nocip_FR_{win_width}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_nocip_FR_{win_width}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6    
    else:
        plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)   
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_div_nocip_FR_{win_width}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_div_nocip_FR_{win_width}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6            
    
    plt.show()
    plt.close()    
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+4-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
       
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Anchor sets below.
    if set_name=='Lr_5.5_FR_FE': 
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_R'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Lr_12.5_FR_FE':    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_R'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Se_5.5_FR_FE':    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_R'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        
    elif set_name=='Se_12.5_FR_FE':    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_R'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_R'][::-1])),       linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07            
    elif set_name=='Lr_5.5_FR_vicinity_FE':     
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_RN'][::-1])) / (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_RN'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR near GCS ({Chi_sets["Cfx_0.01_Lr_5.5_cip_FN"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_RnN'][::-1])),    linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR not near GCS ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_12.5_FR_vicinity_FE':   
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_RN'][::-1])) / (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR near Chi ({Chi_sets["Cfx_0.01_Lr_12.5_cip_FN"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR not near Chi ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_12.5_FR_vicinity_FE':    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_RN'][::-1])) / (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR near Chi ({Chi_sets["Cfx_0.01_Se_12.5_cip_FN"]+Chi_sets["Cfx_0.01_Se_12.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR not near Chi ({Chi_sets["Cfx_0.01_Se_12.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_5.5_FR_vicinity_FE':    
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_RN'][::-1])) / (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_FN']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_RN'][::-1])),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR near GCS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FN"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RN"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_RnN'][::-1]))/ (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_FnN'])+ np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_RnN'][::-1])),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR not near GCS ({Chi_sets["Cfx_0.01_Se_5.5_nocip_FnN"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_RnN"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_5.5_FR_Cl_FE':   
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_R'][::-1]))   /(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Lr_12.5_FR_Cl_FE':   
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_R'][::-1]))   /(np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_5.5_FR_Cl_FE':     
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_F'])  + np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_R'][::-1]))   /(np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_F'])  + np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_12.5_FR_Cl_FE':   
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_RCUS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_FCUS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_RCUS'][::-1])),  linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR Chi close to GCS in US ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCUS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCUS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_RCDS'][::-1]))/(np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_FCDS'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_RCDS'][::-1])),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR Chi close to GCS in DS ({Chi_sets["Cfx_0.01_Se_5.5_cip_FCDS"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RCDS"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
        plot1.plot(positions_sm, (np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_R'][::-1]))   /(np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_F'])   +np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_R'][::-1])),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip FR all Chi ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=1) #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Lr_5.5_US_vs_OTH':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_US']) /np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip US ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_OTH'])/np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_US']) /np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip US ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_OTH'])/np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_US']) /np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip US ({Chi_sets["Cfx_0.01_Se_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_OTH'])/np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Se_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_US']) /np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip US ({Chi_sets["Cfx_0.01_Se_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_OTH'])/np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip/nocip OTH ({Chi_sets["Cfx_0.01_Se_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    
    elif set_name=='Lr_5.5_US_vs_OTH_nocip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip DS ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH_nocip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip DS ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH_nocip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip DS ({Chi_sets["Cfx_0.01_Se_5.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip OTH ({Chi_sets["Cfx_0.01_Se_5.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH_nocip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip DS ({Chi_sets["Cfx_0.01_Se_12.5_nocip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip OTH ({Chi_sets["Cfx_0.01_Se_12.5_nocip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    
    elif set_name=='Lr_5.5_US_vs_OTH_cip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip DS ({Chi_sets["Cfx_0.01_Lr_5.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip OTH ({Chi_sets["Cfx_0.01_Lr_5.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Lr_12.5_US_vs_OTH_cip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip DS ({Chi_sets["Cfx_0.01_Lr_12.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip OTH ({Chi_sets["Cfx_0.01_Lr_12.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    elif set_name=='Se_5.5_US_vs_OTH_cip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 5.5 cip US ({Chi_sets["Cfx_0.01_Se_5.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 5.5 cip OTH ({Chi_sets["Cfx_0.01_Se_5.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07       
    elif set_name=='Se_12.5_US_vs_OTH_cip':      
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_US']),   linestyle='-',  color='#757d8b', linewidth=2.5, alpha=1, label=f'Se 12.5 cip DS ({Chi_sets["Cfx_0.01_Se_12.5_cip_US"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_OTH']),  linestyle='-',  color='#e4d1b4', linewidth=2.5, alpha=1, label=f'Se 12.5 cip OTH ({Chi_sets["Cfx_0.01_Se_12.5_cip_OTH"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07         
    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+8+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    plot1.tick_params(axis='both', labelsize=20)
    #plot1.set_xticks([0, length], minor='True')
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_xlim([-30000, 30008])
    if set_name in ['Se_12.5_US_vs_OTH_cip', 'Se_5.5_US_vs_OTH_cip', 'Lr_12.5_US_vs_OTH_cip', 'Lr_5.5_US_vs_OTH_cip']:
        plot1.set_ylabel(f'RPKM cip', size=20)
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6                   
    elif set_name in ['Se_12.5_US_vs_OTH_nocip', 'Se_5.5_US_vs_OTH_nocip', 'Lr_12.5_US_vs_OTH_nocip', 'Lr_5.5_US_vs_OTH_nocip']:
        plot1.set_ylabel(f'RPKM nocip', size=20)  
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_nocip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_nocip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6                    
    else:
        plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_div_nocip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_cip_div_nocip_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp_lim_30000.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6            
    
    plt.show()
    plt.close()  
   
    return


#Lr Se Cfx effects.

#Output path.
Out_path=f'{PWD}\Chi_anchor_analysis_GCS_vicinity_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Plot_combinations\\'
Dir_check_create(Out_path)
#Name of the signal to plotted (protein or smth.).
Signal_name='Lr_5.5_FR_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_FR_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_FR_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_FR_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_5.5_FR_vicinity_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_FR_vicinity_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_FR_vicinity_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_FR_vicinity_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)

Out_path=f'{PWD}\Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Plot_combinations\\'
Dir_check_create(Out_path)
Signal_name='Lr_5.5_FR_Cl_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_FR_Cl_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_FR_Cl_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_FR_Cl_FE'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)

Out_path=f'{PWD}\Chi_anchor_analysis_US_Chi_vs_others\\50000_width_50_smoothing_1000_vicinity_del_cor\Plot_combinations\\'
Dir_check_create(Out_path)
Signal_name='Lr_5.5_US_vs_OTH'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_US_vs_OTH'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_US_vs_OTH'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_US_vs_OTH'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)

Signal_name='Lr_5.5_US_vs_OTH_cip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_US_vs_OTH_cip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_US_vs_OTH_cip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_US_vs_OTH_cip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)

Signal_name='Lr_5.5_US_vs_OTH_nocip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_US_vs_OTH_nocip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_US_vs_OTH_nocip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_US_vs_OTH_nocip'
plot_FE_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)



#######
#Plot the signal for different groups of anchors together.
#Coverage depth or RPKM.
#######


def plot_RPKM_over_anchors(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes in sets.
    Chi_sets={'Cfx_0.01_Lr_5.5_cip_F'    :    395,
              'Cfx_0.01_Lr_5.5_cip_R'    :    438,
              'Cfx_0.01_Lr_12.5_cip_F'   :    395,
              'Cfx_0.01_Lr_12.5_cip_R'   :    438,
              'Cfx_0.01_Lr_5.5_nocip_F'  :    395,
              'Cfx_0.01_Lr_5.5_nocip_R'  :    438,
              'Cfx_0.01_Lr_12.5_nocip_F' :    395,
              'Cfx_0.01_Lr_12.5_nocip_R' :    438,
              'Cfx_0.01_Se_5.5_cip_F'    :    395,
              'Cfx_0.01_Se_5.5_cip_R'    :    438,
              'Cfx_0.01_Se_12.5_cip_F'   :    395,
              'Cfx_0.01_Se_12.5_cip_R'   :    438,
              'Cfx_0.01_Se_5.5_nocip_F'  :    395,
              'Cfx_0.01_Se_5.5_nocip_R'  :    438,
              'Cfx_0.01_Se_12.5_nocip_F' :    395,
              'Cfx_0.01_Se_12.5_nocip_R' :    438,
              'Cfx_0.01_Lr_5.5_cip_FR'    :   395,
              'Cfx_0.01_Lr_5.5_cip_RF'    :   438,
              'Cfx_0.01_Lr_12.5_cip_FR'   :   395,
              'Cfx_0.01_Lr_12.5_cip_RF'   :   438,
              'Cfx_0.01_Lr_5.5_nocip_FR'  :   395,
              'Cfx_0.01_Lr_5.5_nocip_RF'  :   438,
              'Cfx_0.01_Lr_12.5_nocip_FR' :   395,
              'Cfx_0.01_Lr_12.5_nocip_RF' :   438,
              'Cfx_0.01_Se_5.5_cip_FR'    :   395,
              'Cfx_0.01_Se_5.5_cip_RF'    :   438,
              'Cfx_0.01_Se_12.5_cip_FR'   :   395,
              'Cfx_0.01_Se_12.5_cip_RF'   :   438,
              'Cfx_0.01_Se_5.5_nocip_FR'  :   395,
              'Cfx_0.01_Se_5.5_nocip_RF'  :   438,
              'Cfx_0.01_Se_12.5_nocip_FR' :   395,
              'Cfx_0.01_Se_12.5_nocip_RF' :   438,              
              }
    
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=30000
    length=5000
    for name, file in wig_in_dict.items():
        data=wig_FE_over_genes_parsing(name, file)
        win_width=data[1]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+4, 1)    
    print(positions[0], positions[-1])
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    ##Anchor sets below.
    if set_name=='Lr_5.5_FR':
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_F']) + np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_R'][::-1]),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_R'][::-1]),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_R"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_12.5_FR':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_F']) + np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_R'][::-1]),   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip FR ({Chi_sets["Cfx_0.01_Lr_12.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_R'][::-1]), linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip FR ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_F"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_R"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_5.5_FR':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_F']) + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_R'][::-1]),     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip FR ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_R'][::-1]),   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip FR ({Chi_sets["Cfx_0.01_Se_5.5_nocip_F"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_R"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Se_12.5_FR':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_F']) + np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_R'][::-1]),   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip FR ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_F'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_R'][::-1]), linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip FR ({Chi_sets["Cfx_0.01_Se_12.5_nocip_F"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_R"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
    elif set_name=='Lr_5.5':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_F'])  +  np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_R'][::-1]),      linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip F ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+ Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_FR']) +  np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_cip_RF'][::-1]),     linestyle='--', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip R ({Chi_sets["Cfx_0.01_Lr_5.5_cip_FR"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_F']) + np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_R'][::-1]),    linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip F ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_F"]+ Chi_sets["Cfx_0.01_Lr_5.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_FR'])+ np.array(dict_of_wigs['Cfx_0.01_Lr_5.5_nocip_RF'][::-1]),   linestyle='--', color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip R ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_FR"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Lr_12.5':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_F'])  + np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_R'][::-1]),    linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip F ({Chi_sets["Cfx_0.01_Lr_12.5_cip_F"]+ Chi_sets["Cfx_0.01_Lr_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_FR']) + np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_cip_RF'][::-1]),   linestyle='--', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip R ({Chi_sets["Cfx_0.01_Lr_12.5_cip_FR"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_F']) +np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_R'][::-1]),  linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip F ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_F"]+ Chi_sets["Cfx_0.01_Lr_12.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_FR'])+np.array(dict_of_wigs['Cfx_0.01_Lr_12.5_nocip_RF'][::-1]), linestyle='--', color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip R ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_FR"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_5.5':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_F'])  +  np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_R'][::-1]),      linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip F ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+ Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_FR']) +  np.array(dict_of_wigs['Cfx_0.01_Se_5.5_cip_RF'][::-1]),     linestyle='--', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip R ({Chi_sets["Cfx_0.01_Se_5.5_cip_FR"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_F']) + np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_R'][::-1]),    linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip F ({Chi_sets["Cfx_0.01_Se_5.5_nocip_F"]+ Chi_sets["Cfx_0.01_Se_5.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_FR'])+ np.array(dict_of_wigs['Cfx_0.01_Se_5.5_nocip_RF'][::-1]),   linestyle='--', color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip R ({Chi_sets["Cfx_0.01_Se_5.5_nocip_FR"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
    elif set_name=='Se_12.5':
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_F'])  + np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_R'][::-1]),      linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip F ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+ Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_FR']) + np.array(dict_of_wigs['Cfx_0.01_Se_12.5_cip_RF'][::-1]),     linestyle='--', color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip R ({Chi_sets["Cfx_0.01_Se_12.5_cip_FR"]+Chi_sets["Cfx_0.01_Se_12.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_F']) +np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_R'][::-1]),    linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip F ({Chi_sets["Cfx_0.01_Se_12.5_nocip_F"]+ Chi_sets["Cfx_0.01_Se_12.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.plot(positions, np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_FR'])+np.array(dict_of_wigs['Cfx_0.01_Se_12.5_nocip_RF'][::-1]),   linestyle='--', color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip R ({Chi_sets["Cfx_0.01_Se_12.5_nocip_FR"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
          
                     
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+8+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    plot1.tick_params(axis='both', labelsize=20)
    #plot1.set_xticks([0, length], minor='True')
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM', size=20)
    plot1.set_xlim([-35000, 35008])
    plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_RPKM_FR_{win_width}bp.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_RPKM_FR_{win_width}bp.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6    
    plt.show()
    plt.close()   
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+4-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
       
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Anchor sets below.
    if set_name=='Lr_5.5_FR':
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        Lr_5_5_cip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_R'][::-1])
        Lr_5_5_cip_ar_sc=Lr_5_5_cip_ar/np.mean(np.concatenate([Lr_5_5_cip_ar[:15000], Lr_5_5_cip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Lr_5_5_cip_ar_sc,     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip FR ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})') #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.annotate('', xy=(10000,np.mean(Lr_5_5_cip_ar_sc[60000:70000])), xytext=(10000,np.min(Lr_5_5_cip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_5_5_cip_delta=np.mean(Lr_5_5_cip_ar_sc[60000:70000])-np.min(Lr_5_5_cip_ar_sc)
        Lr_5_5_cip_text=r'$\Delta$='+str(round(Lr_5_5_cip_delta,2))
        plot1.annotate(Lr_5_5_cip_text, xy=(10000,(np.mean(Lr_5_5_cip_ar_sc[60000:70000])+np.min(Lr_5_5_cip_ar_sc))/2), xytext=(10000,(np.mean(Lr_5_5_cip_ar_sc[60000:70000])+np.min(Lr_5_5_cip_ar_sc))/2), zorder=11)
        
        Lr_5_5_nocip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_R'][::-1])
        Lr_5_5_nocip_ar_sc=Lr_5_5_nocip_ar/np.mean(np.concatenate([Lr_5_5_nocip_ar[:15000], Lr_5_5_nocip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Lr_5_5_nocip_ar_sc,   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip FR ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_F"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_R"]})') #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17             
        plot1.annotate('', xy=(15000,np.mean(Lr_5_5_nocip_ar_sc[60000:70000])), xytext=(15000,np.min(Lr_5_5_nocip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=12)
        Lr_5_5_nocip_delta=np.mean(Lr_5_5_nocip_ar_sc[60000:70000])-np.min(Lr_5_5_nocip_ar_sc)
        Lr_5_5_nocip_text=r'$\Delta$='+str(round(Lr_5_5_nocip_delta,2)) 
        plot1.annotate(Lr_5_5_nocip_text, xy=(15000,(np.mean(Lr_5_5_nocip_ar_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_sc))/2), xytext=(15000,(np.mean(Lr_5_5_nocip_ar_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_sc))/2), zorder=13)
    
    elif set_name=='Lr_12.5_FR':   
        Lr_12_5_cip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_F']) + np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_R'][::-1])
        Lr_12_5_cip_ar_sc=Lr_12_5_cip_ar/np.mean(np.concatenate([Lr_12_5_cip_ar[:15000], Lr_12_5_cip_ar[-15000:]], axis=None))     
        plot1.plot(positions_sm, Lr_12_5_cip_ar_sc,   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip FR ({Chi_sets["Cfx_0.01_Lr_12.5_cip_F"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_R"]})') #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.annotate('', xy=(10000,np.mean(Lr_12_5_cip_ar_sc[60000:70000])), xytext=(10000,np.min(Lr_12_5_cip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_cip_delta=np.mean(Lr_12_5_cip_ar_sc[60000:70000])-np.min(Lr_12_5_cip_ar_sc)
        Lr_12_5_cip_text=r'$\Delta$='+str(round(Lr_12_5_cip_delta,2))
        plot1.annotate(Lr_12_5_cip_text, xy=(10000,(np.mean(Lr_12_5_cip_ar_sc[60000:70000])+np.min(Lr_12_5_cip_ar_sc))/2), xytext=(10000,(np.mean(Lr_12_5_cip_ar_sc[60000:70000])+np.min(Lr_12_5_cip_ar_sc))/2), zorder=11)
        
        Lr_12_5_nocip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_R'][::-1])
        Lr_12_5_nocip_ar_sc=Lr_12_5_nocip_ar/np.mean(np.concatenate([Lr_12_5_nocip_ar[:15000], Lr_12_5_nocip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Lr_12_5_nocip_ar_sc, linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip FR ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_F"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_R"]})') #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
        plot1.annotate('', xy=(15000,np.mean(Lr_12_5_nocip_ar_sc[60000:70000])), xytext=(15000,np.min(Lr_12_5_nocip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_nocip_delta=np.mean(Lr_12_5_nocip_ar_sc[60000:70000])-np.min(Lr_12_5_nocip_ar_sc)
        Lr_12_5_nocip_text=r'$\Delta$='+str(round(Lr_12_5_nocip_delta,2))
        plot1.annotate(Lr_12_5_nocip_text, xy=(15000,(np.mean(Lr_12_5_nocip_ar_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_sc))/2), xytext=(15000,(np.mean(Lr_12_5_nocip_ar_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_sc))/2), zorder=11)

    elif set_name=='Se_5.5_FR':
        Se_5_5_cip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_F']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_R'][::-1])
        Se_5_5_cip_ar_sc=Se_5_5_cip_ar/np.mean(np.concatenate([Se_5_5_cip_ar[:15000], Se_5_5_cip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Se_5_5_cip_ar_sc,     linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip FR ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})') #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.annotate('', xy=(10000,np.mean(Se_5_5_cip_ar_sc[60000:70000])), xytext=(10000,np.min(Se_5_5_cip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Se_5_5_cip_delta=np.mean(Se_5_5_cip_ar_sc[60000:70000])-np.min(Se_5_5_cip_ar_sc)
        Se_5_5_cip_text=r'$\Delta$='+str(round(Se_5_5_cip_delta,2))
        plot1.annotate(Se_5_5_cip_text, xy=(10000,(np.mean(Se_5_5_cip_ar_sc[60000:70000])+np.min(Se_5_5_cip_ar_sc))/2), xytext=(10000,(np.mean(Se_5_5_cip_ar_sc[60000:70000])+np.min(Se_5_5_cip_ar_sc))/2))
        
        Se_5_5_nocip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_R'][::-1])
        Se_5_5_nocip_ar_sc=Se_5_5_nocip_ar/np.mean(np.concatenate([Se_5_5_nocip_ar[:15000], Se_5_5_nocip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Se_5_5_nocip_ar_sc,   linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip FR ({Chi_sets["Cfx_0.01_Se_5.5_nocip_F"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_R"]})') #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
        plot1.annotate('', xy=(15000,np.mean(Se_5_5_nocip_ar_sc[60000:70000])), xytext=(15000,np.min(Se_5_5_nocip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=12)
        Se_5_5_nocip_delta=np.mean(Se_5_5_nocip_ar_sc[60000:70000])-np.min(Se_5_5_nocip_ar_sc)
        Se_5_5_nocip_text=r'$\Delta$='+str(round(Se_5_5_nocip_delta,2)) 
        plot1.annotate(Se_5_5_nocip_text, xy=(15000,(np.mean(Se_5_5_nocip_ar_sc[60000:70000])+np.min(Se_5_5_nocip_ar_sc))/2), xytext=(15000,(np.mean(Se_5_5_nocip_ar_sc[60000:70000])+np.min(Se_5_5_nocip_ar_sc))/2), zorder=13)        
    
    elif set_name=='Se_12.5_FR':
        Se_12_5_cip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_F']) + np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_R'][::-1])
        Se_12_5_cip_ar_sc=Se_12_5_cip_ar/np.mean(np.concatenate([Se_12_5_cip_ar[:15000], Se_12_5_cip_ar[-15000:]], axis=None))   
        plot1.plot(positions_sm, Se_12_5_cip_ar_sc,   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip FR ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})') #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.annotate('', xy=(10000,np.mean(Se_12_5_cip_ar_sc[60000:70000])), xytext=(10000,np.min(Se_12_5_cip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Se_12_5_cip_delta=np.mean(Se_12_5_cip_ar_sc[60000:70000])-np.min(Se_12_5_cip_ar_sc)
        Se_12_5_cip_text=r'$\Delta$='+str(round(Se_12_5_cip_delta,2))
        plot1.annotate(Se_12_5_cip_text, xy=(10000,(np.mean(Se_12_5_cip_ar_sc[60000:70000])+np.min(Se_12_5_cip_ar_sc))/2), xytext=(10000,(np.mean(Se_12_5_cip_ar_sc[60000:70000])+np.min(Se_12_5_cip_ar_sc))/2))
        
        Se_12_5_nocip_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_R'][::-1])
        Se_12_5_nocip_ar_sc=Se_12_5_nocip_ar/np.mean(np.concatenate([Se_12_5_nocip_ar[:15000], Se_12_5_nocip_ar[-15000:]], axis=None))
        plot1.plot(positions_sm, Se_12_5_nocip_ar_sc, linestyle='-',  color='#333738', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip FR ({Chi_sets["Cfx_0.01_Se_12.5_nocip_F"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_R"]})') #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17      
        plot1.annotate('', xy=(15000,np.mean(Se_12_5_nocip_ar_sc[60000:70000])), xytext=(15000,np.min(Se_12_5_nocip_ar_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Se_12_5_nocip_delta=np.mean(Se_12_5_nocip_ar_sc[60000:70000])-np.min(Se_12_5_nocip_ar_sc)
        Se_12_5_nocip_text=r'$\Delta$='+str(round(Se_12_5_nocip_delta,2))
        plot1.annotate(Se_12_5_nocip_text, xy=(15000,(np.mean(Se_12_5_nocip_ar_sc[60000:70000])+np.min(Se_12_5_nocip_ar_sc))/2), xytext=(15000,(np.mean(Se_12_5_nocip_ar_sc[60000:70000])+np.min(Se_12_5_nocip_ar_sc))/2))
    
    elif set_name=='Lr_5.5_cip':
        Lr_5_5_cip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_R'][::-1])
        Lr_5_5_cip_ar_F_sc=Lr_5_5_cip_F_ar/np.mean(np.concatenate([Lr_5_5_cip_F_ar[:15000], Lr_5_5_cip_F_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Lr_5_5_cip_F_ar,      linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip F ({Chi_sets["Cfx_0.01_Lr_5.5_cip_F"]+ Chi_sets["Cfx_0.01_Lr_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(8000,np.mean(Lr_5_5_cip_ar_F_sc[60000:70000])), xytext=(8000,np.min(Lr_5_5_cip_ar_F_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_5_5_cip_F_delta=round((np.mean(Lr_5_5_cip_ar_F_sc[60000:70000])-np.min(Lr_5_5_cip_ar_F_sc))*100,0)
        Lr_5_5_cip_F_text=r'$\Delta$='+str(Lr_5_5_cip_F_delta)+'%'
        print(f'Lr 5.5 cip FF and RR: {Lr_5_5_cip_F_text}')
        #plot1.annotate(Lr_5_5_cip_F_text, xy=(8000,(np.mean(Lr_5_5_cip_ar_F_sc[60000:70000])+np.min(Lr_5_5_cip_ar_F_sc))/2), xytext=(8000,(np.mean(Lr_5_5_cip_ar_F_sc[60000:70000])+np.min(Lr_5_5_cip_ar_F_sc))/2), zorder=11)
        
        Lr_5_5_cip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_cip_RF'][::-1])
        Lr_5_5_cip_ar_FR_sc=Lr_5_5_cip_FR_ar/np.mean(np.concatenate([Lr_5_5_cip_FR_ar[:15000], Lr_5_5_cip_FR_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Lr_5_5_cip_FR_ar,     linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Lr 5.5 cip R ({Chi_sets["Cfx_0.01_Lr_5.5_cip_FR"]+Chi_sets["Cfx_0.01_Lr_5.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(12000,np.mean(Lr_5_5_cip_ar_FR_sc[60000:70000])), xytext=(12000,np.min(Lr_5_5_cip_ar_FR_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_5_5_cip_FR_delta=round((np.mean(Lr_5_5_cip_ar_FR_sc[60000:70000])-np.min(Lr_5_5_cip_ar_FR_sc))*100,0)
        Lr_5_5_cip_FR_text=r'$\Delta$='+str(Lr_5_5_cip_FR_delta)+'%'
        print(f'Lr 5.5 cip FR and RF: {Lr_5_5_cip_FR_text}')
        #plot1.annotate(Lr_5_5_cip_FR_text, xy=(12000,(np.mean(Lr_5_5_cip_ar_FR_sc[60000:70000])+np.min(Lr_5_5_cip_ar_FR_sc))/2), xytext=(12000,(np.mean(Lr_5_5_cip_ar_FR_sc[60000:70000])+np.min(Lr_5_5_cip_ar_FR_sc))/2), zorder=11)
    
    elif set_name=='Lr_5.5_nocip':
        Lr_5_5_nocip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_R'][::-1])
        Lr_5_5_nocip_ar_F_sc=Lr_5_5_nocip_F_ar/np.mean(np.concatenate([Lr_5_5_nocip_F_ar[:15000], Lr_5_5_nocip_F_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Lr_5_5_nocip_F_ar,    linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip F ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_F"]+ Chi_sets["Cfx_0.01_Lr_5.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(16000,np.mean(Lr_5_5_nocip_ar_F_sc[60000:70000])), xytext=(16000,np.min(Lr_5_5_nocip_ar_F_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_5_5_nocip_F_delta=round((np.mean(Lr_5_5_nocip_ar_F_sc[60000:70000])-np.min(Lr_5_5_nocip_ar_F_sc))*100,0)
        Lr_5_5_nocip_F_text=r'$\Delta$='+str(Lr_5_5_nocip_F_delta)+'%'
        print(f'Lr 5.5 nocip FF and RR: {Lr_5_5_nocip_F_text}')
        #plot1.annotate(Lr_5_5_nocip_F_text, xy=(16000,(np.mean(Lr_5_5_nocip_ar_F_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_F_sc))/2), xytext=(16000,(np.mean(Lr_5_5_nocip_ar_F_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_F_sc))/2), zorder=11)
        
        Lr_5_5_nocip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_5.5_nocip_RF'][::-1])
        Lr_5_5_nocip_ar_FR_sc=Lr_5_5_nocip_FR_ar/np.mean(np.concatenate([Lr_5_5_nocip_FR_ar[:15000], Lr_5_5_nocip_FR_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Lr_5_5_nocip_FR_ar,   linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Lr 5.5 nocip R ({Chi_sets["Cfx_0.01_Lr_5.5_nocip_FR"]+Chi_sets["Cfx_0.01_Lr_5.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(20000,np.mean(Lr_5_5_nocip_ar_FR_sc[60000:70000])), xytext=(20000,np.min(Lr_5_5_nocip_ar_FR_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_5_5_nocip_FR_delta=round((np.mean(Lr_5_5_nocip_ar_FR_sc[60000:70000])-np.min(Lr_5_5_nocip_ar_FR_sc))*100,0)
        Lr_5_5_nocip_FR_text=r'$\Delta$='+str(Lr_5_5_nocip_FR_delta)+'%'
        print(f'Lr 5.5 nocip FR and RF: {Lr_5_5_nocip_FR_text}')
        #plot1.annotate(Lr_5_5_nocip_FR_text, xy=(20000,(np.mean(Lr_5_5_nocip_ar_FR_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_FR_sc))/2), xytext=(20000,(np.mean(Lr_5_5_nocip_ar_FR_sc[60000:70000])+np.min(Lr_5_5_nocip_ar_FR_sc))/2), zorder=11)
        
    elif set_name=='Lr_12.5_cip':
        Lr_12_5_cip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_R'][::-1])
        Lr_12_5_cip_ar_F_sc=Lr_12_5_cip_F_ar/np.mean(np.concatenate([Lr_12_5_cip_F_ar[:15000], Lr_12_5_cip_F_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Lr_12_5_cip_F_ar,    linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip F ({Chi_sets["Cfx_0.01_Lr_12.5_cip_F"]+ Chi_sets["Cfx_0.01_Lr_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(8000,np.mean(Lr_12_5_cip_ar_F_sc[60000:70000])), xytext=(8000,np.min(Lr_12_5_cip_ar_F_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_cip_F_delta=round((np.mean(Lr_12_5_cip_ar_F_sc[60000:70000])-np.min(Lr_12_5_cip_ar_F_sc))*100,0)
        Lr_12_5_cip_F_text=r'$\Delta$='+str(Lr_12_5_cip_F_delta)+'%'
        print(f'Lr 12.5 cip FF and RR: {Lr_12_5_cip_F_text}')
        #plot1.annotate(Lr_12_5_cip_F_text, xy=(8000,(np.mean(Lr_12_5_cip_ar_F_sc[60000:70000])+np.min(Lr_12_5_cip_ar_F_sc))/2), xytext=(8000,(np.mean(Lr_12_5_cip_ar_F_sc[60000:70000])+np.min(Lr_12_5_cip_ar_F_sc))/2), zorder=11)
        
        Lr_12_5_cip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_cip_RF'][::-1])
        Lr_12_5_cip_ar_FR_sc=Lr_12_5_cip_FR_ar/np.mean(np.concatenate([Lr_12_5_cip_FR_ar[:15000], Lr_12_5_cip_FR_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Lr_12_5_cip_FR_ar,   linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Lr 12.5 cip R ({Chi_sets["Cfx_0.01_Lr_12.5_cip_FR"]+Chi_sets["Cfx_0.01_Lr_12.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(12000,np.mean(Lr_12_5_cip_ar_FR_sc[60000:70000])), xytext=(12000,np.min(Lr_12_5_cip_ar_FR_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_cip_FR_delta=round((np.mean(Lr_12_5_cip_ar_FR_sc[60000:70000])-np.min(Lr_12_5_cip_ar_FR_sc))*100,0)
        Lr_12_5_cip_FR_text=r'$\Delta$='+str(Lr_12_5_cip_FR_delta)+'%'
        print(f'Lr 12.5 cip FR and RF: {Lr_12_5_cip_FR_text}')
        #plot1.annotate(Lr_12_5_cip_FR_text, xy=(12000,(np.mean(Lr_12_5_cip_ar_FR_sc[60000:70000])+np.min(Lr_12_5_cip_ar_FR_sc))/2), xytext=(12000,(np.mean(Lr_12_5_cip_ar_FR_sc[60000:70000])+np.min(Lr_12_5_cip_ar_FR_sc))/2), zorder=11)
   
    elif set_name=='Lr_12.5_nocip':     
        Lr_12_5_nocip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_R'][::-1])
        Lr_12_5_nocip_ar_F_sc=Lr_12_5_nocip_F_ar/np.mean(np.concatenate([Lr_12_5_nocip_F_ar[:15000], Lr_12_5_nocip_F_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Lr_12_5_nocip_F_ar,  linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip F ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_F"]+ Chi_sets["Cfx_0.01_Lr_12.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(16000,np.mean(Lr_12_5_nocip_ar_F_sc[60000:70000])), xytext=(16000,np.min(Lr_12_5_nocip_ar_F_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_nocip_F_delta=round((np.mean(Lr_12_5_nocip_ar_F_sc[60000:70000])-np.min(Lr_12_5_nocip_ar_F_sc))*100,0)
        Lr_12_5_nocip_F_text=r'$\Delta$='+str(Lr_12_5_nocip_F_delta)+'%'
        print(f'Lr 12.5 nocip FF and RR: {Lr_12_5_nocip_F_text}')
        #plot1.annotate(Lr_12_5_nocip_F_text, xy=(16000,(np.mean(Lr_12_5_nocip_ar_F_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_F_sc))/2), xytext=(16000,(np.mean(Lr_12_5_nocip_ar_F_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_F_sc))/2), zorder=11)
        
        Lr_12_5_nocip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Lr_12.5_nocip_RF'][::-1])
        Lr_12_5_nocip_ar_FR_sc=Lr_12_5_nocip_FR_ar/np.mean(np.concatenate([Lr_12_5_nocip_FR_ar[:15000], Lr_12_5_nocip_FR_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Lr_12_5_nocip_FR_ar, linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Lr 12.5 nocip R ({Chi_sets["Cfx_0.01_Lr_12.5_nocip_FR"]+Chi_sets["Cfx_0.01_Lr_12.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.annotate('', xy=(20000,np.mean(Lr_12_5_nocip_ar_FR_sc[60000:70000])), xytext=(20000,np.min(Lr_12_5_nocip_ar_FR_sc)), arrowprops=dict(arrowstyle='<->', facecolor='black'), zorder=10)
        Lr_12_5_nocip_FR_delta=round((np.mean(Lr_12_5_nocip_ar_FR_sc[60000:70000])-np.min(Lr_12_5_nocip_ar_FR_sc))*100,0)
        Lr_12_5_nocip_FR_text=r'$\Delta$='+str(Lr_12_5_nocip_FR_delta)+'%'
        print(f'Lr 12.5 nocip FR and RF: {Lr_12_5_nocip_FR_text}')
        #plot1.annotate(Lr_12_5_nocip_FR_text, xy=(20000,(np.mean(Lr_12_5_nocip_ar_FR_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_FR_sc))/2), xytext=(20000,(np.mean(Lr_12_5_nocip_ar_FR_sc[60000:70000])+np.min(Lr_12_5_nocip_ar_FR_sc))/2), zorder=11)
        
    elif set_name=='Se_5.5_cip':
        Se_5_5_cip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_R'][::-1])
        Se_5_5_cip_ar_F_sc=Se_5_5_cip_F_ar/np.mean(np.concatenate([Se_5_5_cip_F_ar[:15000], Se_5_5_cip_F_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Se_5_5_cip_F_ar,      linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Se 5.5 cip F ({Chi_sets["Cfx_0.01_Se_5.5_cip_F"]+ Chi_sets["Cfx_0.01_Se_5.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_5_5_cip_F_delta=round((np.mean(Se_5_5_cip_ar_F_sc[60000:70000])-np.min(Se_5_5_cip_ar_F_sc))*100,0)
        Se_5_5_cip_F_text=r'$\Delta$='+str(Se_5_5_cip_F_delta)+'%'
        print(f'Se 5.5 cip FF and RR: {Se_5_5_cip_F_text}')
        
        Se_5_5_cip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_cip_RF'][::-1])
        Se_5_5_cip_ar_FR_sc=Se_5_5_cip_FR_ar/np.mean(np.concatenate([Se_5_5_cip_FR_ar[:15000], Se_5_5_cip_FR_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Se_5_5_cip_FR_ar,     linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Se 5.5 cip R ({Chi_sets["Cfx_0.01_Se_5.5_cip_FR"]+Chi_sets["Cfx_0.01_Se_5.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_5_5_cip_FR_delta=round((np.mean(Se_5_5_cip_ar_FR_sc[60000:70000])-np.min(Se_5_5_cip_ar_FR_sc))*100,0)
        Se_5_5_cip_FR_text=r'$\Delta$='+str(Se_5_5_cip_FR_delta)+'%'
        print(f'Se 5.5 cip FR and RF: {Se_5_5_cip_FR_text}')
        
    elif set_name=='Se_5.5_nocip':   
        Se_5_5_nocip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_R'][::-1])
        Se_5_5_nocip_ar_F_sc=Se_5_5_nocip_F_ar/np.mean(np.concatenate([Se_5_5_nocip_F_ar[:15000], Se_5_5_nocip_F_ar[-15000:]], axis=None))       
        plot1.plot(positions_sm, Se_5_5_nocip_F_ar,    linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip F ({Chi_sets["Cfx_0.01_Se_5.5_nocip_F"]+ Chi_sets["Cfx_0.01_Se_5.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_5_5_nocip_F_delta=round((np.mean(Se_5_5_nocip_ar_F_sc[60000:70000])-np.min(Se_5_5_nocip_ar_F_sc))*100,0)
        Se_5_5_nocip_F_text=r'$\Delta$='+str(Se_5_5_nocip_F_delta)+'%'
        print(f'Se 5.5 nocip FF and RR: {Se_5_5_nocip_F_text}')
        
        Se_5_5_nocip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_5.5_nocip_RF'][::-1])
        Se_5_5_nocip_ar_FR_sc=Se_5_5_nocip_FR_ar/np.mean(np.concatenate([Se_5_5_nocip_FR_ar[:15000], Se_5_5_nocip_FR_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Se_5_5_nocip_FR_ar,   linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Se 5.5 nocip R ({Chi_sets["Cfx_0.01_Se_5.5_nocip_FR"]+Chi_sets["Cfx_0.01_Se_5.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_5_5_nocip_FR_delta=round((np.mean(Se_5_5_nocip_ar_FR_sc[60000:70000])-np.min(Se_5_5_nocip_ar_FR_sc))*100,0)
        Se_5_5_nocip_FR_text=r'$\Delta$='+str(Se_5_5_nocip_FR_delta)+'%'
        print(f'Se 5.5 nocip FR and RF: {Se_5_5_nocip_FR_text}')
        
    elif set_name=='Se_12.5_cip':
        Se_12_5_cip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_R'][::-1])
        Se_12_5_cip_ar_F_sc=Se_12_5_cip_F_ar/np.mean(np.concatenate([Se_12_5_cip_F_ar[:15000], Se_12_5_cip_F_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Se_12_5_cip_F_ar,      linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Se 12.5 cip F ({Chi_sets["Cfx_0.01_Se_12.5_cip_F"]+ Chi_sets["Cfx_0.01_Se_12.5_cip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_12_5_cip_F_delta=round((np.mean(Se_12_5_cip_ar_F_sc[60000:70000])-np.min(Se_12_5_cip_ar_F_sc))*100,0)
        Se_12_5_cip_F_text=r'$\Delta$='+str(Se_12_5_cip_F_delta)+'%'
        print(f'Se 12.5 cip FF and RR: {Se_12_5_cip_F_text}')
        
        Se_12_5_cip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_cip_RF'][::-1])
        Se_12_5_cip_ar_FR_sc=Se_12_5_cip_FR_ar/np.mean(np.concatenate([Se_12_5_cip_FR_ar[:15000], Se_12_5_cip_FR_ar[-15000:]], axis=None))        
        plot1.plot(positions_sm, Se_12_5_cip_FR_ar,     linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Se 12.5 cip R ({Chi_sets["Cfx_0.01_Se_12.5_cip_FR"]+Chi_sets["Cfx_0.01_Se_12.5_cip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_12_5_cip_FR_delta=round((np.mean(Se_12_5_cip_ar_FR_sc[60000:70000])-np.min(Se_12_5_cip_ar_FR_sc))*100,0)
        Se_12_5_cip_FR_text=r'$\Delta$='+str(Se_12_5_cip_FR_delta)+'%'
        print(f'Se 12.5 cip FR and RF: {Se_12_5_cip_FR_text}')
        
    elif set_name=='Se_12.5_nocip':
        Se_12_5_nocip_F_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_F'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_R'][::-1])
        Se_12_5_nocip_ar_F_sc=Se_12_5_nocip_F_ar/np.mean(np.concatenate([Se_12_5_nocip_F_ar[:15000], Se_12_5_nocip_F_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Se_12_5_nocip_F_ar,    linestyle='-',  color='#42634A', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip F ({Chi_sets["Cfx_0.01_Se_12.5_nocip_F"]+ Chi_sets["Cfx_0.01_Se_12.5_nocip_R"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        Se_12_5_nocip_F_delta=round((np.mean(Se_12_5_nocip_ar_F_sc[60000:70000])-np.min(Se_12_5_nocip_ar_F_sc))*100,0)
        Se_12_5_nocip_F_text=r'$\Delta$='+str(Se_12_5_nocip_F_delta)+'%'
        print(f'Se 12.5 nocip FF and RR: {Se_12_5_nocip_F_text}')
        
        Se_12_5_nocip_FR_ar=np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_FR'])+np.array(dict_of_wigs_sm['Cfx_0.01_Se_12.5_nocip_RF'][::-1])
        Se_12_5_nocip_ar_FR_sc=Se_12_5_nocip_FR_ar/np.mean(np.concatenate([Se_12_5_nocip_FR_ar[:15000], Se_12_5_nocip_FR_ar[-15000:]], axis=None))         
        plot1.plot(positions_sm, Se_12_5_nocip_FR_ar,   linestyle='-', color='#BDBDBD', linewidth=2.5, alpha=1, label=f'Se 12.5 nocip R ({Chi_sets["Cfx_0.01_Se_12.5_nocip_FR"]+Chi_sets["Cfx_0.01_Se_12.5_nocip_RF"]})', zorder=7) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07        Se_12_5_nocip_FR_delta=round((np.mean(Se_12_5_nocip_ar_FR_sc[60000:70000])-np.min(Se_12_5_nocip_ar_FR_sc))*100,0)
        Se_12_5_nocip_FR_text=r'$\Delta$='+str(Se_12_5_nocip_FR_delta)+'%'
        print(f'Se 12.5 nocip FR and RF: {Se_12_5_nocip_FR_text}')        
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+8+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    plot1.tick_params(axis='both', labelsize=20)
    #plot1.set_xticks([0, length], minor='True')
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM', size=20)
    plot1.set_xlim([-35000, 35008])
    plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_RPKM_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.savefig(f'{output_path}\\{set_name}_over_{set_type}_RPKM_FR_smoothed_{win_width}bp_smoothed_{2*sm_window}bp.svg', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6    
    plt.show()
    plt.close()  
   
    return


#Lr Se Cfx effects, investigate signal on F and R strands.
Out_path=f'{PWD}\Chi_anchor_analysis_Closest_Chi_polyshed\\50000_width_50_smoothing_1000_vicinity_del_cor\Plot_combinations\\'
Dir_check_create(Out_path)
Signal_name='Lr_5.5_FR'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_FR'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_FR'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_FR'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_5.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_5.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Lr_12.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_5.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_cip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)
Signal_name='Se_12.5_nocip'
plot_RPKM_over_anchors(Wig_data_in_dict, Sm_window, Out_path, Signal_name, Set_type)