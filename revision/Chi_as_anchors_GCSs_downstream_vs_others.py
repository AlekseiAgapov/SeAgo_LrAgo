###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#The script identifies closest Chi sites to the GCSs to the 3' of them,
#and compares their enrichment to enrichment at other Chi sites.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import re
import scipy
from scipy import stats
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
from operator import itemgetter

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\"

#Chi-site sequence.
Chi={'F' : 'GCTGGTGG', 'R' : 'CCACCAGC'}

#Path to GCSs files dict.
GCSs_dict={'Cfx_0.01' : 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Cfx_0.01_trusted_GCSs.txt', 
           }

#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\"

#Path to Ago-Seq files.

Dict_of_wigs_path={'Lr_5.5_cip'    : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_plus_FPKM_av.wig',    'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_minus_FPKM_av.wig'},
                   'Lr_5.5_nocip'  : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_plus_FPKM_av.wig',  'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_minus_FPKM_av.wig'},
                   'Lr_12.5_cip'   : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_cip_plus_FPKM_av.wig',   'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_cip_minus_FPKM_av.wig'},
                   'Lr_12.5_nocip' : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_nocip_plus_FPKM_av.wig', 'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_12.5_nocip_minus_FPKM_av.wig'},
                   'Se_5.5_cip'    : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_cip_plus_FPKM_av.wig',    'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_cip_minus_FPKM_av.wig'},
                   'Se_5.5_nocip'  : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_nocip_plus_FPKM_av.wig',  'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_5.5_nocip_minus_FPKM_av.wig'},
                   'Se_12.5_cip'   : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_cip_plus_FPKM_av.wig',   'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_cip_minus_FPKM_av.wig'},
                   'Se_12.5_nocip' : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_nocip_plus_FPKM_av.wig', 'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Se_12.5_nocip_minus_FPKM_av.wig'},                                                                               
                   }

#Dict_of_wigs_path={'Lr_5.5_cip'    : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_plus_FPKM_av.wig',    'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_cip_minus_FPKM_av.wig'},
#                   'Lr_5.5_nocip'  : {'F': PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_plus_FPKM_av.wig',  'R' : PWD + 'Ago_Seq\Syn_Lro_libraries_data_wig_RPKM_polyshed_av\Lr_5.5_nocip_minus_FPKM_av.wig'},
#                   }

#Path to file with deletions.
Deletions_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Deletions_NC_012971_2.broadPeak"

#Path to regions to be masked (like rRNA operons).
Masking_regs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\Regions_to_mask_polysh.bed"

#Path to the genome FASTA.
Genome_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\Scripts\Additional_genome_features\\E_coli_BL21_DE3.fasta"

#Distance threshold between GCS and chi site, bp.
Vicinity_thr=1000

#Region width.
Region_width=50000

#Smoothin window width.
Smoothing_win_width=50

#Output path.
Output_path=f'{PWD}Chi_anchor_analysis_US_Chi_vs_others\\{Region_width}_width_{Smoothing_win_width}_smoothing_{Vicinity_thr}_vicinity_del_cor\\'

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
Dir_check_create(Output_path+'Signal_of_GCSs_wig_full\\') 


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
#Genome sequence parsing.
#######

def genome_seq(genome_path):
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genome_sequence=str(record.seq)
    genome.close()
    print('Whole genome average GC: ' + str(SeqUtils.GC(genome_sequence)))
    print('Whole genome length: ' + str(len(genome_sequence)))        
    return genome_sequence


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
#Get coordinates of Chi sites.
#######

def get_chi_coords(Chi_seq_dict, deletions, genome_sequence):
    
    Chi_for_seq=Chi_seq_dict['F']
    Chi_rev_seq=Chi_seq_dict['R']

    Chi_for_coords=[m.start() for m in re.finditer(Chi_for_seq, genome_sequence)]
    Chi_for_coords_del_cor=[]
    for chi_site in Chi_for_coords:
        delited=0
        for deletion in deletions:
            if deletion[1]>=chi_site>=deletion[0]:
                delited=1
        if delited==0:            
            Chi_for_coords_del_cor.append(chi_site)    
            
    Chi_rev_coords=[m.start() for m in re.finditer(Chi_rev_seq, genome_sequence)]
    Chi_rev_coords_del_cor=[]
    for chi_site in Chi_rev_coords:
        delited=0
        for deletion in deletions:
            if deletion[1]>=chi_site>=deletion[0]:
                delited=1
        if delited==0:            
            Chi_rev_coords_del_cor.append(chi_site)      
    
    print(f'Number of Chi sites of positive strand: {len(Chi_for_coords_del_cor)}')
    print(f'Number of Chi sites of negative strand: {len(Chi_rev_coords_del_cor)}')
    
    return Chi_for_coords_del_cor, Chi_rev_coords_del_cor


#######
#Get coordinates of Chi sites near GCSs.
#Get coordinates of Chi sites without GCSs nearby.
#######

def Chi_sites_and_GCSs(chi_coords_ar, chi_orientation, GCSs_coords_ar, vicinity_thr):
        
    chi_closest_to_GCS_us=[]
    chi_closest_to_GCS_ds=[] 
    
    GCSs_info_dict={}
    Chis_info_dict={}
    
    for GCS_site in GCSs_coords_ar:
        
        vicinity_us=10000000
        vicinity_ds=10000000 
        
        for chi_site in chi_coords_ar:
            #US: a Chi-site is upstream of a GCS (GCS is downstream of Chi); DS: a Chi-site is downstream of a GCS (GCS is upstream of Chi).
            if chi_orientation=='plus':
                if 0<=GCS_site-chi_site<=vicinity_us:
                    vicinity_us=GCS_site-chi_site
                    chi_site_us=chi_site
                if 0<=chi_site-GCS_site<=vicinity_ds:
                    vicinity_ds=chi_site-GCS_site
                    chi_site_ds=chi_site
            if chi_orientation=='minus':
                if 0<=GCS_site-chi_site<=vicinity_ds:
                    vicinity_ds=GCS_site-chi_site
                    chi_site_ds=chi_site
                if 0<=chi_site-GCS_site<=vicinity_us:
                    vicinity_us=chi_site-GCS_site
                    chi_site_us=chi_site
        
        if ('chi_site_us' in locals()) and ('chi_site_ds' in locals()):  
            chi_closest_to_GCS_us.append(chi_site_us)
            chi_closest_to_GCS_ds.append(chi_site_ds)
            
            GCSs_info_dict[GCS_site]={'US' : [chi_site_us, vicinity_us, chi_orientation], 'DS' : [chi_site_ds, vicinity_ds, chi_orientation]}
            
            if chi_site_us not in Chis_info_dict:
                Chis_info_dict[chi_site_us]={'US': [[GCS_site, vicinity_us, chi_orientation]]}
            else:
                if 'US' in Chis_info_dict[chi_site_us]:
                    Chis_info_dict[chi_site_us]['US'].append([GCS_site, vicinity_us, chi_orientation])
                else:
                    Chis_info_dict[chi_site_us]['US']=[[GCS_site, vicinity_us, chi_orientation]]
            
            if chi_site_ds not in Chis_info_dict:
                Chis_info_dict[chi_site_ds]={'DS': [[GCS_site, vicinity_ds, chi_orientation]]}
            else:
                if 'DS' in Chis_info_dict[chi_site_ds]:
                    Chis_info_dict[chi_site_ds]['DS'].append([GCS_site, vicinity_ds, chi_orientation]) 
                else:
                    Chis_info_dict[chi_site_ds]['DS']=[[GCS_site, vicinity_ds, chi_orientation]]
            
            print(f'For GCS position: {GCS_site} and Chi orientation "{chi_orientation}", closest Chi US position: {chi_site_us}; closest Chi DS position: {chi_site_ds}')
    
    #Get non-redundant sets of Chi-sites.
    chi_closest_to_GCS_us_unique=list(set(chi_closest_to_GCS_us))
    print(f'Total number of Chi sites {len(chi_closest_to_GCS_us)} in US; number of unique sites {len(chi_closest_to_GCS_us_unique)} in US')
    chi_closest_to_GCS_ds_unique=list(set(chi_closest_to_GCS_ds))
    print(f'Total number of Chi sites {len(chi_closest_to_GCS_ds)} in DS; number of unique sites {len(chi_closest_to_GCS_ds_unique)} in DS')
    print(chi_closest_to_GCS_us_unique)
    print(chi_closest_to_GCS_ds_unique)
    chi_closest_to_GCS_usds_intersec=[x for x in chi_closest_to_GCS_us_unique if x in chi_closest_to_GCS_ds_unique]
    print(f'Number of Chi sites both in US and DS: {len(chi_closest_to_GCS_usds_intersec)}, which is {float(len(chi_closest_to_GCS_usds_intersec))*100/len(chi_closest_to_GCS_us_unique)}% of US and {float(len(chi_closest_to_GCS_usds_intersec))*100/len(chi_closest_to_GCS_ds_unique)}% of DS')
    
    #Get unabmbigously classfied Chi-sites.
    chi_closest_to_GCS_us_unique_only=[x for x in chi_closest_to_GCS_us_unique if x not in chi_closest_to_GCS_usds_intersec]
    chi_closest_to_GCS_ds_unique_only=[x for x in chi_closest_to_GCS_ds_unique if x not in chi_closest_to_GCS_usds_intersec]
    print(f'Number of unique and filtered Chi sites in US {len(chi_closest_to_GCS_us_unique_only)};  Number of unique and filtered Chi sites in DS {len(chi_closest_to_GCS_ds_unique_only)}\n')
    
    #Get non-us Chi-sites (other than US).
    Not_US_chi_ar=[]
    for chi_site in chi_coords_ar:
        if chi_site not in chi_closest_to_GCS_us_unique_only:
            Not_US_chi_ar.append(chi_site)
 
    return chi_closest_to_GCS_us_unique_only, Not_US_chi_ar


#######
#Collect Ago signal under Chi-sites.
#######

def collect_Ago_signal(Chi_list, strand, Ago_wig, win_width):
    
    Chi_RPKM_data_tab={}
    Chi_RPKM_data_wig={}
    Chi_metasignal=np.array([0.0]*(win_width+win_width+4))
    for Chi_coord in Chi_list:
        if (Chi_coord-win_width>0) and (Chi_coord+win_width+4)<len(Ago_wig):
            if strand=='plus':
                Signal_around_Chi=np.array(Ago_wig[Chi_coord-win_width:Chi_coord+win_width+4])
                Chi_metasignal+=Signal_around_Chi
                Chi_RPKM_data_tab[Chi_coord]=[np.mean(Ago_wig[Chi_coord-win_width:Chi_coord]), np.mean(Ago_wig[Chi_coord:Chi_coord+win_width+4])]
                Chi_RPKM_data_wig[Chi_coord]=np.array(Ago_wig[Chi_coord-5000:Chi_coord+5000+4])
            elif strand=='minus':
                Signal_around_Chi=np.array(Ago_wig[Chi_coord-win_width:Chi_coord+win_width+4][::-1])
                Chi_metasignal+=Signal_around_Chi
                Chi_RPKM_data_tab[Chi_coord]=[np.mean(Ago_wig[Chi_coord-win_width:Chi_coord][::-1]), np.mean(Ago_wig[Chi_coord:Chi_coord+win_width+4][::-1])]
                Chi_RPKM_data_wig[Chi_coord]=np.array(Ago_wig[Chi_coord-5000:Chi_coord+5000+4][::-1])                    
    
    return Chi_metasignal, Chi_RPKM_data_tab, Chi_RPKM_data_wig


#######
#Write .tab file with PRKM info for US and DS.
#######

def join_F_and_R_data(Chi_metasignal_F, Chi_metasignal_R, Chi_RPKM_data_F_tab, Chi_RPKM_data_R_tab, Chi_RPKM_data_F_wig, Chi_RPKM_data_R_wig):
    
    Chi_metasignal=(Chi_metasignal_F+Chi_metasignal_R)/(len(Chi_RPKM_data_F_tab)+len(Chi_RPKM_data_R_tab))
    
    Chi_RPKM_data_tab_FR_dict=Chi_RPKM_data_F_tab
    Chi_RPKM_data_tab_R_dict=Chi_RPKM_data_R_tab
    Chi_RPKM_data_tab_FR_dict.update(Chi_RPKM_data_tab_R_dict)
    
    Chi_RPKM_data_wig_FR_dict=Chi_RPKM_data_F_wig
    Chi_RPKM_data_wig_R_dict=Chi_RPKM_data_R_wig
    Chi_RPKM_data_wig_FR_dict.update(Chi_RPKM_data_wig_R_dict)
        
    return Chi_metasignal, Chi_RPKM_data_tab_FR_dict, Chi_RPKM_data_wig_FR_dict


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
#Write fasta-like file.
#######

def write_anchors_full_data(dict_of_tracks, fileout_path):
    fileout=open(fileout_path, 'w')
    for track_name, track_data in dict_of_tracks.items():
        fileout.write(f'{track_name}\t')
        for position in track_data:
            fileout.write(f'{round(position, 2)}, ')
        fileout.write('\n')
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
    plt.close()     
    
    return


#######
#Calculate mean FE, FE STD, and return FE array for further statistics.
#######

def calc_FE_STD(data_cip, data_nocip, set_name):
      
    if len(data_cip)==len(data_nocip):
        print(f'Number of Chi sites considered: {len(data_nocip)}')
    else:
        print(f'Number of Chi sites in cip sample: {len(data_cip)}\nNumber of Chi sites in nocip sample: {len(data_nocip)}\n')
        return
    
    FE_values=[]
    
    for chi_site, ago_data in data_cip.items():
        Site_mean_cip=np.mean(data_cip[chi_site][5000+4:])
        Site_mean_nocip=np.mean(data_nocip[chi_site][5000+4:])
        if Site_mean_nocip>0:
            Site_FE=Site_mean_cip/Site_mean_nocip
        else:
            Site_FE=(Site_mean_cip+1)/(Site_mean_nocip+1)
        FE_values.append(Site_FE)
    
    FE_mean_value=np.mean(FE_values)
    FE_STD_value=np.std(FE_values)

    print(f'For sample {set_name}, FE mean={round(FE_mean_value, 2)}, FE STD={round(FE_STD_value, 2)}')
    
    return FE_mean_value, FE_STD_value, FE_values


#######
#Prepare data for barplots, and make statistics.
#######

def prep_barplot_data_and_stat(Detailed_metasignal_dict, Dataset_name, GCSs_set_name):
    
    Chi_wig_ar_US_cip = Detailed_metasignal_dict[Dataset_name+'_cip'][GCSs_set_name]['US']
    Chi_wig_ar_OTH_cip = Detailed_metasignal_dict[Dataset_name+'_cip'][GCSs_set_name]['OTH']
    Chi_wig_ar_US_nocip = Detailed_metasignal_dict[Dataset_name+'_nocip'][GCSs_set_name]['US']
    Chi_wig_ar_OTH_nocip = Detailed_metasignal_dict[Dataset_name+'_nocip'][GCSs_set_name]['OTH']    
    
    US_FE_mean_value, US_FE_STD_value, US_FE_values = calc_FE_STD(Chi_wig_ar_US_cip, Chi_wig_ar_US_nocip, Dataset_name)
    OTH_FE_mean_value, OTH_FE_STD_value, OTH_FE_values = calc_FE_STD(Chi_wig_ar_OTH_cip, Chi_wig_ar_OTH_nocip, Dataset_name)
    
    print(f'Working with dataset: {Dataset_name}')
    print('Comparing of Chi-sites with a GCS to the 3\' vs other Chi-sites.')
    Intervals_stat=stats.ttest_ind(US_FE_values, OTH_FE_values, equal_var=True)
    print(f'Sample size 1: {len(US_FE_values)}, Sample size 2: {len(OTH_FE_values)}')
    print(f'\nT-test FE, Mean1={round(np.mean(US_FE_values),3)}; Mean2={round(np.mean(OTH_FE_values),3)}\np-value={Intervals_stat[1]}\nt-statistic={Intervals_stat[0]}\n') 

    return US_FE_mean_value, US_FE_STD_value, OTH_FE_mean_value, OTH_FE_STD_value


def plot_barplots(Mean_data_ar, STD_data_ar, output_path, Dataset_names_ar):
    
    #Plot FE over anchors. Chi and distance from GCSs.
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)

    Loc_max_ar=[]
    for i in range(len(Mean_data_ar)):
        Loc_max=Mean_data_ar[i]+STD_data_ar[i]
        Loc_max_ar.append(Loc_max)
    Glob_max=np.max(Loc_max_ar)    
    
    #US: a Chi-site is upstream of a GCS (GCS is downstream of Chi); OTH: other Chi-sites.
    Conditions=[]
    color_list=[]
    X_coords=[]
    for i in range(len(Dataset_names_ar)):
        if i%2==0:
            bar_name=f'DS\n{Dataset_names_ar[i]}'
            Conditions.append(bar_name)
            color_list.append('#c96458')
            X_coords.append(i+1)     
        else:
            bar_name=f'OTH\n{Dataset_names_ar[i]}'
            Conditions.append(bar_name)
            color_list.append('#6a65c7')
            X_coords.append(i+1)
    
    X_coords_m=X_coords
      
    Bars=plot_av.bar(X_coords, Mean_data_ar, yerr=STD_data_ar, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.9, color=color_list, edgecolor='k', linewidth=1)
    plot_av.set_ylabel('RPKM cip/RPKM nocip', size=16)
    plot_av.set_xticks(X_coords_m)
    plot_av.set_xticklabels(Conditions, rotation=0, size=8)  
    plot_av.set_ylim([0.47, Glob_max*1.1])
    plt.tight_layout()
    plt.savefig(f'{output_path}\\Barplot_DS_vs_others_over_Chi_cip_div_nocip_FR_5000_bp_US.png', dpi=400, figsize=(5, 3), transparent=True)   
    plt.savefig(f'{output_path}\\Barplot_DS_vs_others_over_Chi_cip_div_nocip_FR_5000_bp_US.svg', dpi=400, figsize=(5, 3), transparent=True)   
    plt.close() 
    
    return

#######
#Wrapper function.
#######

def wrapper_function(Chi_seq_dict, dict_of_wigs_path, GCSs_path_dict, deletions_path, masking_regs_path, vicinity_thr, win_width, sm_window, genome_path, output_path):
    
    #Read deletions and regions to be masked.
    deletions=deletions_info(deletions_path)
    masking_regs=deletions_info(masking_regs_path)   
    
    #Reads input data in wig files.
    dict_of_wigs={}
    for Ago_dataset_name, pair_path in dict_of_wigs_path.items():
        wig_data_F=wig_parsing(pair_path['F'], masking_regs)
        wig_data_R=wig_parsing(pair_path['R'], masking_regs)
        dict_of_wigs[Ago_dataset_name]={'F' : wig_data_F, 'R' : wig_data_R}     
    
    #Read genome file, get chi-sites coordinates.
    genome_sequence=genome_seq(genome_path)
    Chi_for_coords, Chi_rev_coords=get_chi_coords(Chi_seq_dict, deletions, genome_sequence)
    
    #Reads gcss input data.
    dict_of_gcss_coords={}
    for gcss_name, gcss_file in GCSs_path_dict.items():
        gcss_data=Read_GCSs(gcss_file, deletions)
        dict_of_gcss_coords[gcss_name]=gcss_data   
        
    #Find closest pairs Chi-site:GCS.
    #Chi_coords_dict={'All' : {'F' : Chi_for_coords, 'R' : Chi_rev_coords}}
    Chi_coords_dict={}
    for gcss_name, gcss_coords in dict_of_gcss_coords.items(): 
        
        Chi_for_closest_to_GCS_us_ar, Chi_for_others_ar=Chi_sites_and_GCSs(Chi_for_coords, 'plus', gcss_coords, vicinity_thr)
        Chi_rev_closest_to_GCS_us_ar, Chi_rev_others_ar=Chi_sites_and_GCSs(Chi_rev_coords, 'minus', gcss_coords, vicinity_thr)
        
        Chi_coords_dict[gcss_name]={'F_CUS' : Chi_for_closest_to_GCS_us_ar, 'F_OTH' : Chi_for_others_ar,
                                    'R_CUS' : Chi_rev_closest_to_GCS_us_ar, 'R_OTH' : Chi_rev_others_ar,}            
    
    Cumulative_metasignal_dict={}
    Detailed_metasignal_dict={}
    for Ago_dataset_name, Ago_pair_wig in dict_of_wigs.items():
        Ago_wig_F=Ago_pair_wig['F']
        Ago_wig_R=Ago_pair_wig['R']
        
        Cumulative_metasignal_dict[Ago_dataset_name]={}
        Detailed_metasignal_dict[Ago_dataset_name]={}        
        
        for GCSs_set_name, Chi_coords_lists in Chi_coords_dict.items():
            print(f'Now processing: {Ago_dataset_name} and {GCSs_set_name}')
    
            Chi_coords_list_FCUS=Chi_coords_lists['F_CUS']
            Chi_coords_list_FOTH=Chi_coords_lists['F_OTH']
            Chi_coords_list_RCUS=Chi_coords_lists['R_CUS']                
            Chi_coords_list_ROTH=Chi_coords_lists['R_OTH']
            
            #Get metasignal.
            Chi_metasignal_FCUS, Chi_RPKM_data_FCUS, Chi_RPKM_data_FCUS_wig_ar=collect_Ago_signal(Chi_coords_list_FCUS, 'plus', Ago_wig_F, win_width)
            Chi_metasignal_FOTH, Chi_RPKM_data_FOTH, Chi_RPKM_data_FOTH_wig_ar=collect_Ago_signal(Chi_coords_list_FOTH, 'plus', Ago_wig_F, win_width)
            Chi_metasignal_RCUS, Chi_RPKM_data_RCUS, Chi_RPKM_data_RCUS_wig_ar=collect_Ago_signal(Chi_coords_list_RCUS, 'minus', Ago_wig_R, win_width)
            Chi_metasignal_ROTH, Chi_RPKM_data_ROTH, Chi_RPKM_data_ROTH_wig_ar=collect_Ago_signal(Chi_coords_list_ROTH, 'minus', Ago_wig_R, win_width)
    
            #Join F and R strands' data.
            Chi_metasignal_US,  Chi_RPKM_data_tab_US,  Chi_RPKM_data_wig_US  = join_F_and_R_data(Chi_metasignal_FCUS, Chi_metasignal_RCUS, Chi_RPKM_data_FCUS, Chi_RPKM_data_RCUS, Chi_RPKM_data_FCUS_wig_ar, Chi_RPKM_data_RCUS_wig_ar)
            Chi_metasignal_OTH, Chi_RPKM_data_tab_OTH, Chi_RPKM_data_wig_OTH = join_F_and_R_data(Chi_metasignal_FOTH, Chi_metasignal_ROTH, Chi_RPKM_data_FOTH, Chi_RPKM_data_ROTH, Chi_RPKM_data_FOTH_wig_ar, Chi_RPKM_data_ROTH_wig_ar)
    
            #Store data.
            Cumulative_metasignal_dict[Ago_dataset_name][GCSs_set_name]={'US' : Chi_metasignal_US, 'OTH' : Chi_metasignal_OTH}
            Detailed_metasignal_dict[Ago_dataset_name][GCSs_set_name] = {'US' : Chi_RPKM_data_wig_US, 'OTH' : Chi_RPKM_data_wig_OTH}
            
            #Write data.
            anchor_plot(Chi_metasignal_US, win_width, sm_window, Ago_dataset_name, GCSs_set_name+'_US', output_path)
            print(f'Writing RPKM over US Chi-sites...')
            write_wig(Chi_metasignal_US, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_US_width_{win_width}bp.wig', f'{win_width}')
            print(f'Writing detailed averaged RPKM for US Chi-sites...')
            write_GSCs_RPKM(Chi_RPKM_data_tab_US, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_US_{win_width}bp.txt')
            print(f'Writing detailed RPKM over US Chi-sites, for each Chi...')
            write_anchors_full_data(Chi_RPKM_data_wig_US, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_US_width_5000bp.txt')                
            
            anchor_plot(Chi_metasignal_OTH, win_width, sm_window, Ago_dataset_name, GCSs_set_name+'_OTH', output_path)
            print(f'Writing RPKM over OTHER Chi-sites...')
            write_wig(Chi_metasignal_OTH, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_OTH_width_{win_width}bp.wig', f'{win_width}')
            print(f'Writing detailed averaged RPKM for OTHER Chi-sites...')
            write_GSCs_RPKM(Chi_RPKM_data_tab_OTH, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_OTH_{win_width}bp.txt')
            print(f'Writing detailed RPKM over OTHER Chi-sites, for each Chi...')
            write_anchors_full_data(Chi_RPKM_data_wig_OTH, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_over_{GCSs_set_name}_OTH_width_5000bp.txt')                     
                         
    
    ##Make barplots.
    ##Anchor sets below.
    
    Dataset_name='Lr_5.5'
    GCSs_set_name='Cfx_0.01'  
    US_FE_mean_value_1, US_FE_STD_value_1, OTH_FE_mean_value_1, OTH_FE_STD_value_1 = prep_barplot_data_and_stat(Detailed_metasignal_dict, Dataset_name, GCSs_set_name)
        
    Dataset_name='Lr_12.5'
    GCSs_set_name='Cfx_0.01'
    US_FE_mean_value_2, US_FE_STD_value_2, OTH_FE_mean_value_2, OTH_FE_STD_value_2 = prep_barplot_data_and_stat(Detailed_metasignal_dict, Dataset_name, GCSs_set_name)
     
    Dataset_name='Se_5.5'
    GCSs_set_name='Cfx_0.01'
    US_FE_mean_value_3, US_FE_STD_value_3, OTH_FE_mean_value_3, OTH_FE_STD_value_3 = prep_barplot_data_and_stat(Detailed_metasignal_dict, Dataset_name, GCSs_set_name)
      
    Dataset_name='Se_12.5'
    GCSs_set_name='Cfx_0.01'
    US_FE_mean_value_4, US_FE_STD_value_4, OTH_FE_mean_value_4, OTH_FE_STD_value_4 = prep_barplot_data_and_stat(Detailed_metasignal_dict, Dataset_name, GCSs_set_name)
    
    Means_ar=[US_FE_mean_value_3, OTH_FE_mean_value_3, US_FE_mean_value_4, OTH_FE_mean_value_4, US_FE_mean_value_1, OTH_FE_mean_value_1, US_FE_mean_value_2, OTH_FE_mean_value_2]
    STD_ar=[US_FE_STD_value_3, OTH_FE_STD_value_3, US_FE_STD_value_4, OTH_FE_STD_value_4, US_FE_STD_value_1, OTH_FE_STD_value_1, US_FE_STD_value_2, OTH_FE_STD_value_2]
    Dataset_names=['Se_5.5', 'Se_5.5', 'Se_12.5', 'Se_12.5', 'Lr_5.5', 'Lr_5.5', 'Lr_12.5', 'Lr_12.5']
    
    plot_barplots(Means_ar, STD_ar, output_path, Dataset_names)   
        
    return

wrapper_function(Chi, Dict_of_wigs_path, GCSs_dict, Deletions_path, Masking_regs_path, Vicinity_thr, Region_width, Smoothing_win_width, Genome_path, Output_path)
