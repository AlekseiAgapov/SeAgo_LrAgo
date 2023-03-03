###############################################
##Dmitry Sutormin, 2022##
##Ago-Seq analysis##

#The script makes an anchor plot for Ago enrichment using GCSs positions as anchors.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
from scipy.stats import binom_test
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils import GC as GC_count
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import AlignIO, motifs
from Bio.Align import AlignInfo
import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\"

#Path to GCSs files dict.
GCSs_dict_1={'Cfx_0.01' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Cfx_0.01_trusted_GCSs.txt', 
             }
GCSs_dict_2={'Cfx_0.01_top' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Cfx_0.01_trusted_GCSs_top_200_N3E.txt', 
             'Cfx_0.01_bot' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Cfx_0.01_trusted_GCSs_bottom_200_N3E.txt',
             }

#Path to Ago-Seq files.
Dict_of_wigs_path={'Lr_5.5_cip_F'   : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_plus_FPKM_av.wig',    
                   'Lr_5.5_cip_R'   : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_minus_FPKM_av.wig',
                   'Lr_5.5_nocip_F' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_plus_FPKM_av.wig',  
                   'Lr_5.5_nocip_R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_minus_FPKM_av.wig',
                   'Lr_12.5_cip_F'  : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_cip_plus_FPKM_av.wig',   
                   'Lr_12.5_cip_R'  : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_cip_minus_FPKM_av.wig',
                   'Lr_12.5_nocip_F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_nocip_plus_FPKM_av.wig', 
                   'Lr_12.5_nocip_R': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_nocip_minus_FPKM_av.wig',
                   'Se_5.5_cip_F'   : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_cip_plus_FPKM_av.wig',    
                   'Se_5.5_cip_R'   : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_cip_minus_FPKM_av.wig',
                   'Se_5.5_nocip_F' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_nocip_plus_FPKM_av.wig',  
                   'Se_5.5_nocip_R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_nocip_minus_FPKM_av.wig',
                   'Se_12.5_cip_F'  : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_cip_plus_FPKM_av.wig',   
                   'Se_12.5_cip_R'  : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_cip_minus_FPKM_av.wig',
                   'Se_12.5_nocip_F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_nocip_plus_FPKM_av.wig', 
                   'Se_12.5_nocip_R': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_nocip_minus_FPKM_av.wig',                                                                               
                   }

#Path to file with deletions.
Deletions_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Deletions_NC_012971_2.broadPeak"

#Path to regions to be masked (like rRNA operons).
Masking_regs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Regions_to_mask_polysh.bed"

#Region width.
Region_width=15000

#Smoothin window width.
Smoothing_win_width=200

#Output path.
Output_path=f'{PWD}Anchor_analysis_polyshed\\GCSs_as_anchors_{Region_width}_width_{Smoothing_win_width}_smoothing\\'

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

Dir_check_create(Output_path)
Dir_check_create(Output_path+'Signal_of_GCSs_tab\\')
Dir_check_create(Output_path+'Signal_of_GCSs_wig\\')  

#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile, masking_regs):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    Total_NE=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
            Total_NE+=float(line[0])
    print('Total number of ends: ' + str(Total_NE))
    wigin.close()
    
    Mean_RPKM=np.mean(NE_values)
    for masking_reg in masking_regs:
        masking_reg_len=masking_reg[1]-masking_reg[0]+1
        masking_reg_mask=[Mean_RPKM]*masking_reg_len
        NE_values[masking_reg[0]:masking_reg[1]]=masking_reg_mask
        
    return NE_values

#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar


#######
#Read GCSs files, get GCSs coordinates.
#######

def Read_GCSs(GCSsfile, deletions):
    print('Now is processing: ' + str(GCSsfile))
    GCSsin=open(GCSsfile, 'r')
    GCSs_coordinates=[]
    for line in GCSsin:
        line=line.rstrip().split('\t')
        if line[0] not in ['GCSs_coordinate']:
            GCS_coord=int(line[0])
            GCS_coord_end=GCS_coord+4
            delited=0
            for deletion in deletions:
                if deletion[1]>=GCS_coord>=deletion[0] or deletion[1]>=GCS_coord_end>=deletion[0]:
                    delited=1
            if delited==0:            
                GCSs_coordinates.append(GCS_coord)
    print('Total number of GSCs: ' + str(len(GCSs_coordinates)))
    GCSsin.close()    
    return GCSs_coordinates


#######
#Collect Ago signal under GCSs.
#######

def collect_Ago_signal(GCSs_list, Ago_wig, win_width):
    
    GSCs_RPKM_data={}
    GCSs_metasignal=np.array([0.0]*(win_width+win_width+4))
    for GCSs_coord in GCSs_list:
        if (GCSs_coord-win_width>0) and (GCSs_coord+win_width+4)<len(Ago_wig):
            Signal_around_GCS=np.array(Ago_wig[GCSs_coord-win_width:GCSs_coord+win_width+4])
            GCSs_metasignal+=Signal_around_GCS
            GSCs_RPKM_data[GCSs_coord]=[np.mean(Ago_wig[GCSs_coord-win_width:GCSs_coord]), np.mean(Ago_wig[GCSs_coord:GCSs_coord+win_width+4])]
        
    GCSs_metasignal=GCSs_metasignal/len(GCSs_list)
    
    return GCSs_metasignal, GSCs_RPKM_data


#######
#Write .tab file with PRKM info for US and DS.
#######

def write_GSCs_RPKM(dict_GCSs, RPKM_track_name, anchor_set_name, path_out):
    fileout=open(path_out, 'w')
    fileout.write(f'GCSs_coordinate\t{RPKM_track_name}_{anchor_set_name}_RPKM_US\t{RPKM_track_name}_{anchor_set_name}_RPKM_DS\n')
    for k, v in dict_GCSs.items():
        fileout.write(f'{k}\t{v[0]}\t{v[1]}\n')
    fileout.close()
    return


#######
#Write .wig file.
#######

def write_wig(ar, fileout_path, name):
    fileout=open(fileout_path, 'w')
    fileout.write(f'track type=wiggle_0 name="{name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=NC_007779.1_w3110_Mu start=1 step=1\n')
    for point in ar:
        fileout.write(f'{point}\n')
    fileout.close()
    return


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
#Plot anchor plot.
#######

def anchor_plot(GCSs_metasignal, win_width, sm_window, Ago_dataset_name, GCSs_set_name, output_path):
    
    length=1000
    positions=np.arange(-win_width, win_width+4, 1)
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)   
    plot1.plot(positions, GCSs_metasignal,   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'{Ago_dataset_name} over GCSs {GCSs_set_name}', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
    ticks=np.arange(-win_width,win_width+4+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)    
    plt.savefig(f'{output_path}\\{Ago_dataset_name}_FE_over_{GCSs_set_name}_{win_width}bp_nd_with_body.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close() 
    
    
    #Smoothing.
    GCSs_metasignal_sm=Smoothing(GCSs_metasignal, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+4-(sm_window-1), 1)  
    print(positions_sm[0], positions_sm[-1])
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    plot1.plot(positions_sm, GCSs_metasignal_sm,   linestyle='-',  color='#B6B8BD', linewidth=2.5, alpha=1, label=f'{Ago_dataset_name} over GCSs {GCSs_set_name}', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
    ticks=np.arange(-win_width,win_width+4+1,length).tolist()   
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='GCSs'
    plot1.set_xticks([0, length], minor='True')
    plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'RPKM cip/RPKM nocip', size=20)    
    plt.savefig(f'{output_path}\\{Ago_dataset_name}_FE_over_{GCSs_set_name}_smoothed_{win_width}bp_nd_with_body_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    return


#######
#Wrapper function.
#######

def wrapper_function(gcss_dict, dict_of_wigs_path, deletions_path, masking_regs_path, win_width, sm_window, output_path):
    
    deletions=deletions_info(deletions_path)
    masking_regs=deletions_info(masking_regs_path)
    
    #Reads input data in wig files.
    dict_of_wigs={}
    for file_name, file_path in dict_of_wigs_path.items():
        wig_data=wig_parsing(file_path, masking_regs)
        dict_of_wigs[file_name]=wig_data  
        
    #Reads gcss input data.
    dict_of_gcss_coords={}
    for gcss_name, gcss_file in gcss_dict.items():
        gcss_data=Read_GCSs(gcss_file, deletions)
        dict_of_gcss_coords[gcss_name]=gcss_data
    
    for Ago_dataset_name, Ago_wig in dict_of_wigs.items():
        for GCSs_set_name, GCSs_list in dict_of_gcss_coords.items():
            print(f'Now processing: {Ago_dataset_name} and {GCSs_set_name}')
            GCSs_metasignal, GSCs_RPKM_data=collect_Ago_signal(GCSs_list, Ago_wig, win_width)
            anchor_plot(GCSs_metasignal, win_width, sm_window, Ago_dataset_name, GCSs_set_name, output_path)
            print(f'Writing RPKM over GCSs...')
            write_wig(GCSs_metasignal, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_width_{win_width}bp.wig', f'{win_width}')
            print(f'Writing RPKM for GCSs\' TU, DS...')
            write_GSCs_RPKM(GSCs_RPKM_data, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_{win_width}bp.txt')
            
    return

wrapper_function(GCSs_dict_1, Dict_of_wigs_path, Deletions_path, Masking_regs_path, Region_width, Smoothing_win_width, Output_path)
wrapper_function(GCSs_dict_2, Dict_of_wigs_path, Deletions_path, Masking_regs_path, Region_width, Smoothing_win_width, Output_path)










