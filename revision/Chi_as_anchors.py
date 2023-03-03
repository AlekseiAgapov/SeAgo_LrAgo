###############################################
##Dmitry Sutormin, 2022##
##Ago-Seq analysis##

#The script makes an anchor plot for Ago enrichment using Chi-sites positions as anchors.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import re
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

#Regime of Chi sites selection. "GCS_vicinity" or "Closest_Chi".
Regime="GCS_vicinity"

#Output path.
Output_path=f'{PWD}Chi_anchor_analysis_{Regime}_polyshed\\{Region_width}_width_{Smoothing_win_width}_smoothing_{Vicinity_thr}_vicinity_del_cor\\'

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

def Chi_sites_and_GCSs(chi_coords_ar, chi_orientation, GCSs_coords_ar, vicinity_thr, regime):
    
    if regime=="GCS_vicinity":
    
        chi_near_GCS=[]
        chi_not_near_GCS=[]
        chi_near_GCS_us=[]
        chi_near_GCS_ds=[]
        
        for chi_site in chi_coords_ar:
            vicinity=0
            vicinity_us=0
            vicinity_ds=0
            
            for GCS_site in GCSs_coords_ar:
                if chi_orientation=='plus':
                    if 0<=GCS_site-chi_site<=vicinity_thr:
                        vicinity_us=1
                    if 0>=GCS_site-chi_site>=-vicinity_thr:
                        vicinity_ds=1
                if chi_orientation=='minus':
                    if 0>=GCS_site-chi_site>=-vicinity_thr:
                        vicinity_us=1
                    if 0<=GCS_site-chi_site<=vicinity_thr:
                        vicinity_ds=1
                if np.abs(GCS_site-chi_site)<=vicinity_thr:
                    vicinity=1
            if vicinity==0:
                chi_not_near_GCS.append(chi_site) 
            else:
                chi_near_GCS.append(chi_site)
            if vicinity_us!=0:
                chi_near_GCS_us.append(chi_site)
            if vicinity_ds!=0:
                chi_near_GCS_ds.append(chi_site)
                
        print(f'Number of chi-sites in the vicinity of {vicinity_thr}bp from GCS is: {len(chi_near_GCS)} sites')
        print(f'Number of chi-sites in the vicinity of {vicinity_thr}bp US from GCS is: {len(chi_near_GCS)} sites')
        print(f'Number of chi-sites in the vicinity of {vicinity_thr}bp DS from GCS is: {len(chi_near_GCS)} sites')
        print(f'Number of chi-sites NOT in the vicinity of {vicinity_thr}bp from GCS is: {len(chi_not_near_GCS)} sites')
        
        return chi_near_GCS, chi_not_near_GCS, chi_near_GCS_us, chi_near_GCS_ds
        
    elif regime=="Closest_Chi":
        
        chi_closest_to_GCS_us=[]
        chi_closest_to_GCS_ds=[]        
        
        for GCS_site in GCSs_coords_ar:
            
            vicinity_us=10000000
            vicinity_ds=10000000 
            
            for chi_site in chi_coords_ar:
                
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
                
                print(f'For GCS position: {GCS_site} and Chi orientation "{chi_orientation}", closest Chi US position: {chi_site_us}; closest Chi DS position: {chi_site_ds}')
        chi_closest_to_GCS_us_unique=list(set(chi_closest_to_GCS_us))
        print(f'Total number of Chi sites {len(chi_closest_to_GCS_us)} in US; number of unique sites {len(chi_closest_to_GCS_us_unique)} in US')
        chi_closest_to_GCS_ds_unique=list(set(chi_closest_to_GCS_ds))
        print(f'Total number of Chi sites {len(chi_closest_to_GCS_ds)} in DS; number of unique sites {len(chi_closest_to_GCS_ds_unique)} in DS')
        print(chi_closest_to_GCS_us_unique)
        print(chi_closest_to_GCS_ds_unique)
        chi_closest_to_GCS_usds_intersec=[x for x in chi_closest_to_GCS_us_unique if x in chi_closest_to_GCS_ds_unique]
        print(f'Number of Chi sites both in US and DS: {len(chi_closest_to_GCS_usds_intersec)}, which is {float(len(chi_closest_to_GCS_usds_intersec))*100/len(chi_closest_to_GCS_us_unique)}% of US and {float(len(chi_closest_to_GCS_usds_intersec))*100/len(chi_closest_to_GCS_ds_unique)}% of DS')
        
        chi_closest_to_GCS_us_unique_only=[x for x in chi_closest_to_GCS_us_unique if x not in chi_closest_to_GCS_usds_intersec]
        chi_closest_to_GCS_ds_unique_only=[x for x in chi_closest_to_GCS_ds_unique if x not in chi_closest_to_GCS_usds_intersec]
        print(f'Number of unique and filtered Chi sites in US {len(chi_closest_to_GCS_us_unique_only)};  Number of unique and filtered Chi sites in DS {len(chi_closest_to_GCS_ds_unique_only)}\n')
                          
        return chi_closest_to_GCS_us_unique_only, chi_closest_to_GCS_ds_unique_only


#######
#Collect Ago signal under GCSs.
#######

def collect_Ago_signal(GCSs_list, Ago_wig, win_width):
    
    GSCs_RPKM_data_tab={}
    GSCs_RPKM_data_wig={}
    GCSs_metasignal=np.array([0.0]*(win_width+win_width+4))
    for GCSs_coord in GCSs_list:
        if (GCSs_coord-win_width>0) and (GCSs_coord+win_width+4)<len(Ago_wig):
            Signal_around_GCS=np.array(Ago_wig[GCSs_coord-win_width:GCSs_coord+win_width+4])
            GCSs_metasignal+=Signal_around_GCS
            GSCs_RPKM_data_tab[GCSs_coord]=[np.mean(Ago_wig[GCSs_coord-win_width:GCSs_coord]), np.mean(Ago_wig[GCSs_coord:GCSs_coord+win_width+4])]
            GSCs_RPKM_data_wig[GCSs_coord]=np.array(Ago_wig[GCSs_coord-5000:GCSs_coord+5000+4])
        
    GCSs_metasignal=GCSs_metasignal/len(GCSs_list)
    
    return GCSs_metasignal, GSCs_RPKM_data_tab, GSCs_RPKM_data_wig


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
#Wrapper function.
#######

def wrapper_function(Chi_seq_dict, dict_of_wigs_path, GCSs_path_dict, deletions_path, masking_regs_path, regime, vicinity_thr, win_width, sm_window, genome_path, output_path):
    
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
        
    #Find Chi-site near GCS.
    Chi_coords_dict={'All' : {'F' : Chi_for_coords, 'R' : Chi_rev_coords}}
    for gcss_name, gcss_coords in dict_of_gcss_coords.items():
        if regime=="GCS_vicinity":
            Chi_for_near_GCS, Chi_for_not_near_GCS, Chi_for_near_GCS_us, Chi_for_near_GCS_ds=Chi_sites_and_GCSs(Chi_for_coords, 'plus', gcss_coords, vicinity_thr, regime)
            Chi_rev_near_GCS, Chi_rev_not_near_GCS, Chi_rev_near_GCS_us, Chi_rev_near_GCS_ds=Chi_sites_and_GCSs(Chi_rev_coords, 'minus', gcss_coords, vicinity_thr, regime)
            Chi_coords_dict[gcss_name]={'F_N' : Chi_for_near_GCS, 'F_nN' : Chi_for_not_near_GCS, 'F_NUS' : Chi_for_near_GCS_us, 'F_NDS' : Chi_for_near_GCS_ds,
                                        'R_N' : Chi_rev_near_GCS, 'R_nN' : Chi_rev_not_near_GCS, 'R_NUS' : Chi_rev_near_GCS_us, 'R_NDS' : Chi_rev_near_GCS_ds,}
            
        elif regime=="Closest_Chi":
            Chi_for_closest_to_GCS_us, Chi_for_closest_to_GCS_ds=Chi_sites_and_GCSs(Chi_for_coords, 'plus', gcss_coords, vicinity_thr, regime)
            Chi_rev_closest_to_GCS_us, Chi_rev_closest_to_GCS_ds=Chi_sites_and_GCSs(Chi_rev_coords, 'minus', gcss_coords, vicinity_thr, regime)
            Chi_coords_dict[gcss_name]={'F_CUS' : Chi_for_closest_to_GCS_us, 'F_CDS' : Chi_for_closest_to_GCS_ds,
                                        'R_CUS' : Chi_rev_closest_to_GCS_us, 'R_CDS' : Chi_rev_closest_to_GCS_ds,}            
    
    for Ago_dataset_name, Ago_pair_wig in dict_of_wigs.items():
        Ago_wig_F=Ago_pair_wig['F']
        Ago_wig_R=Ago_pair_wig['R']
        
        for GCSs_set_name, Chi_coords_lists in Chi_coords_dict.items():
            print(f'Now processing: {Ago_dataset_name} and {GCSs_set_name}')
            if GCSs_set_name=='All':      
                Chi_coords_list_F=Chi_coords_lists['F']
                Chi_coords_list_R=Chi_coords_lists['R']
                
                GCSs_metasignal_FF, GSCs_RPKM_data_FF, GSCs_RPKM_data_FF_wig=collect_Ago_signal(Chi_coords_list_F, Ago_wig_F, win_width)
                GCSs_metasignal_RR, GSCs_RPKM_data_RR, GSCs_RPKM_data_RR_wig=collect_Ago_signal(Chi_coords_list_R, Ago_wig_R, win_width)
                GCSs_metasignal_FR, GSCs_RPKM_data_FR, GSCs_RPKM_data_FR_wig=collect_Ago_signal(Chi_coords_list_F, Ago_wig_R, win_width)
                GCSs_metasignal_RF, GSCs_RPKM_data_RF, GSCs_RPKM_data_RF_wig=collect_Ago_signal(Chi_coords_list_R, Ago_wig_F, win_width)                
                
                anchor_plot(GCSs_metasignal_FF, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_F', output_path)
                print(f'Writing RPKM over GCSs...')
                write_wig(GCSs_metasignal_FF, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_F_width_{win_width}bp.wig', f'{win_width}')
                print(f'Writing RPKM for GCSs\' TU, DS...')
                write_GSCs_RPKM(GSCs_RPKM_data_FF, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_F_{win_width}bp.txt')
                print(f'Writing RPKM over GCSs, for each GCS...')
                write_anchors_full_data(GSCs_RPKM_data_FF_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_F_width_5000bp.txt')
                
                anchor_plot(GCSs_metasignal_RR, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_R', output_path)
                print(f'Writing RPKM over GCSs...')
                write_wig(GCSs_metasignal_RR, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_R_width_{win_width}bp.wig', f'{win_width}')
                print(f'Writing RPKM for GCSs\' TU, DS...')
                write_GSCs_RPKM(GSCs_RPKM_data_RR, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_R_{win_width}bp.txt')
                print(f'Writing RPKM over GCSs, for each GCS...')
                write_anchors_full_data(GSCs_RPKM_data_RR_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_R_width_5000bp.txt')                
            
            
                anchor_plot(GCSs_metasignal_FR, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_R', output_path)
                print(f'Writing RPKM over GCSs...')
                write_wig(GCSs_metasignal_FR, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_R_width_{win_width}bp.wig', f'{win_width}')
                print(f'Writing RPKM for GCSs\' TU, DS...')
                write_GSCs_RPKM(GSCs_RPKM_data_FR, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_R_{win_width}bp.txt')
                print(f'Writing RPKM over GCSs, for each GCS...')
                write_anchors_full_data(GSCs_RPKM_data_FR_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_R_width_5000bp.txt')
                
                anchor_plot(GCSs_metasignal_RF, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_F', output_path)
                print(f'Writing RPKM over GCSs...')
                write_wig(GCSs_metasignal_RF, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_F_width_{win_width}bp.wig', f'{win_width}')
                print(f'Writing RPKM for GCSs\' TU, DS...')
                write_GSCs_RPKM(GSCs_RPKM_data_RF, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_F_{win_width}bp.txt')
                print(f'Writing RPKM over GCSs, for each GCS...')
                write_anchors_full_data(GSCs_RPKM_data_RF_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_F_width_5000bp.txt')                
                
            else:
                if regime=="Closest_Chi":
                    
                    Chi_coords_list_FCUS=Chi_coords_lists['F_CUS']
                    Chi_coords_list_FCDS=Chi_coords_lists['F_CDS']
                    Chi_coords_list_RCUS=Chi_coords_lists['R_CUS']                
                    Chi_coords_list_RCDS=Chi_coords_lists['R_CDS']
                    
                    GCSs_metasignal_FCUS, GSCs_RPKM_data_FCUS, GSCs_RPKM_data_FCUS_wig=collect_Ago_signal(Chi_coords_list_FCUS, Ago_wig_F, win_width)
                    GCSs_metasignal_FCDS, GSCs_RPKM_data_FCDS, GSCs_RPKM_data_FCDS_wig=collect_Ago_signal(Chi_coords_list_FCDS, Ago_wig_F, win_width)
                    GCSs_metasignal_RCUS, GSCs_RPKM_data_RCUS, GSCs_RPKM_data_RCUS_wig=collect_Ago_signal(Chi_coords_list_RCUS, Ago_wig_R, win_width)
                    GCSs_metasignal_RCDS, GSCs_RPKM_data_RCDS, GSCs_RPKM_data_RCDS_wig=collect_Ago_signal(Chi_coords_list_RCDS, Ago_wig_R, win_width)
                    
                    anchor_plot(GCSs_metasignal_FCUS, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_FCUS', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_FCUS, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCUS_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_FCUS, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCUS_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_FCUS_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCUS_width_5000bp.txt')                
                    
                    anchor_plot(GCSs_metasignal_FCDS, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_FCDS', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_FCDS, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCDS_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_FCDS, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCDS_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_FCDS_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FCDS_width_5000bp.txt')                     
                    
                    anchor_plot(GCSs_metasignal_RCUS, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_RCUS', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_RCUS, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCUS_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_RCUS, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCUS_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_RCUS_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCUS_width_5000bp.txt')                     
                    
                    anchor_plot(GCSs_metasignal_RCDS, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_RCDS', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_RCDS, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCDS_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_RCDS, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCDS_{win_width}bp.txt')                
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_RCDS_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RCDS_width_5000bp.txt')                     
                           
                elif regime=="GCS_vicinity":
                    
                    Chi_coords_list_FN=Chi_coords_lists['F_N']               
                    Chi_coords_list_FnN=Chi_coords_lists['F_nN']
                    Chi_coords_list_RN=Chi_coords_lists['R_N']                
                    Chi_coords_list_RnN=Chi_coords_lists['R_nN']
                    
                    GCSs_metasignal_FNF,   GSCs_RPKM_data_FNF,   GSCs_RPKM_data_FNF_wig=collect_Ago_signal(Chi_coords_list_FN, Ago_wig_F, win_width)
                    GCSs_metasignal_FnNF,  GSCs_RPKM_data_FnNF,  GSCs_RPKM_data_FnNF_wig=collect_Ago_signal(Chi_coords_list_FnN, Ago_wig_F, win_width)
                    GCSs_metasignal_RNR,   GSCs_RPKM_data_RNR,   GSCs_RPKM_data_RNR_wig=collect_Ago_signal(Chi_coords_list_RN, Ago_wig_R, win_width)
                    GCSs_metasignal_RnNR,  GSCs_RPKM_data_RnNR,  GSCs_RPKM_data_RnNR_wig=collect_Ago_signal(Chi_coords_list_RnN, Ago_wig_R, win_width)
                    
                    anchor_plot(GCSs_metasignal_FNF, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_FN', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_FNF, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FN_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_FNF, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FN_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_FNF_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FN_width_5000bp.txt')                     
                    
                    anchor_plot(GCSs_metasignal_FnNF, win_width, sm_window, Ago_dataset_name+'_F', GCSs_set_name+'_FnN', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_FnNF, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FnN_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_FnNF, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FnN_{win_width}bp.txt')                
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_FnNF_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_F_over_{GCSs_set_name}_FnN_width_5000bp.txt') 
                    
                    anchor_plot(GCSs_metasignal_RNR, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_RN', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_RNR, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RN_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_RNR, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RN_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_RNR_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RN_width_5000bp.txt') 
                    
                    anchor_plot(GCSs_metasignal_RnNR, win_width, sm_window, Ago_dataset_name+'_R', GCSs_set_name+'_RnN', output_path)
                    print(f'Writing RPKM over GCSs...')
                    write_wig(GCSs_metasignal_RnNR, f'{output_path}\Signal_of_GCSs_wig\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RnN_width_{win_width}bp.wig', f'{win_width}')
                    print(f'Writing RPKM for GCSs\' TU, DS...')
                    write_GSCs_RPKM(GSCs_RPKM_data_RnNR, Ago_dataset_name, GCSs_set_name, f'{output_path}\Signal_of_GCSs_tab\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RnN_{win_width}bp.txt')
                    print(f'Writing RPKM over GCSs, for each GCS...')
                    write_anchors_full_data(GSCs_RPKM_data_RnNR_wig, f'{output_path}\Signal_of_GCSs_wig_full\\Signal_{Ago_dataset_name}_R_over_{GCSs_set_name}_RnN_width_5000bp.txt')
                    
    return

wrapper_function(Chi, Dict_of_wigs_path, GCSs_dict, Deletions_path, Masking_regs_path, Regime, Vicinity_thr, Region_width, Smoothing_win_width, Genome_path, Output_path)
