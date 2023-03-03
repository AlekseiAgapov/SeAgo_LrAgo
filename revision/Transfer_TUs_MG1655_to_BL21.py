###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

#Transfer annotation of MG1655 genome with TUs to BL21 genome using blast.
###############################################

#######
#Packages to be imported.
#######

import os
from os import listdir
import numpy as np
import scipy
import pandas as pd
from pandas import DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from matplotlib import cm as cm
from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline


#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\RNA_Seq_BL21\\"

#RegulonDB UTR data
Regulon_UTR="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Databases\\RegulonDB_E_coli\\UTR_5_3_sequence_processed.txt"

#MG1655 U00096.3 sequence
MG1655="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_Gyrase_Topo-Seq\Genomes\\E_coli_K12_MG1655_U00096.3.fasta"

#BL21 genome.
BL21="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Gyrase_cleavage_activity_and_Agos\Genome\\E_coli_BL21_DE3.fasta"

#Path to TUs sequences.
UTRs_sequences_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\RNA_Seq_BL21\\UTRs_full_sequences.fasta"

#Blast results.
Blast_results="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\RNA_Seq_BL21\\UTRs_full_sequences_BL21.blastResult"

#Deletions.
Deletions="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\RNA_Seq_BL21\\Deletions_NC_012971_2.broadPeak"


#######
#Checks if directory exists and if not - creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return


#######
#Deletions parser (BED).
#######

def Deletions_parser(Deletions_inpath):
    Deletions_in=open(Deletions_inpath, 'r')
    Deletions_coord=[]
    for line in Deletions_in:
        line=line.rstrip().split('\t')
        Deletions_coord.append([int(line[1])+1, int(line[2])]) #[TAD start, TAD end]
    Deletions_in.close()
    return Deletions_coord

Deletions_data=Deletions_parser(Deletions)

#######
#Read table with UTR data.
#######

def read_UTR(utr_path):
    #UTR data is +1-based (according to RegulonDB description).
    UTRs=pd.read_csv(utr_path, sep='\t')
    UTRs_coordinates_list=UTRs['UTR_coordinates'].tolist()
    TSS=[]
    TES=[]
    for UTR in UTRs_coordinates_list:
        UTR=UTR.split('-')
        TSS.append(UTR[0])
        TES.append(UTR[1])
    UTRs['TSS']=TSS
    UTRs['TES']=TES
    #print(UTRs)           
    return UTRs

UTRs_info=read_UTR(Regulon_UTR)


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
#Return TUs sequences.
#######

def return_TU_seqs(MG1655_genome_path, UTRs_seqs_path, UTRs):
    MG_genome=genome_seq(MG1655_genome_path)
    UTR_seq_file=open(UTRs_seqs_path, 'w')
    TSS=UTRs['TSS'].tolist()
    TES=UTRs['TES'].tolist()
    Strands=UTRs['Strand'].tolist()
    Names=UTRs['TU_name'].tolist()
    UTRs_seqs=[]
    UTRs_dict={}
    for i in range(len(TSS)):
        seq=MG_genome[int(TSS[i])-1:int(TES[i])]
        UTRs_seqs.append(seq)
        UTR_seq_file.write(f'>{int(TSS[i])}__{TES[i]}__{Strands[i]}__{Names[i]}\n{seq}\n')
        UTRs_dict[f'{int(TSS[i])}__{TES[i]}__{Strands[i]}__{Names[i]}']=seq
    #print(UTRs_seqs[2])
    UTRs['UTR_sequence']=UTRs_seqs
    return UTRs, UTRs_dict, UTRs_seqs

UTRs_info_seq, UTRs_info_seq_dict, UTRs_sequences=return_TU_seqs(MG1655, UTRs_sequences_path, UTRs_info)

#######
#Blast of UTRs performd via command line, because python-implemented blast is too slow.
#blastn -query C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\\RNA_Seq_BL21\\UTRs_full_sequences.fasta -subject C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Gyrase_cleavage_activity_and_Agos\Genome\\E_coli_BL21_DE3.fasta -outfmt=6 -out C:\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\RNA_Seq_BL21\UTRs_full_sequences_BL21.blastResult -max_target_seqs 1 -max_hsps 1
#######

def parse_blast_results(Blast_results_path, UTRs_dict):
    blast_file=open(Blast_results_path, 'r')
    Reversal=-1
    UTRs_BL21_dict={}
    for line in blast_file:
        line=line.rstrip().split('\t')
        query_name=line[0]
        if line[8]!='NaN':
            if int(line[8])<int(line[9]):
                UTR_BL21_start=int(line[8])
                UTR_BL21_end=int(line[9])
                Reversal=0
            elif int(line[8])>int(line[9]):
                UTR_BL21_start=int(line[9])
                UTR_BL21_end=int(line[8])   
                Reversal=1   
        
        UTRs_BL21_dict[query_name]=[UTR_BL21_start, UTR_BL21_end, Reversal]
    
    print(len(UTRs_BL21_dict)) 
    print(len(UTRs_dict))
               
    for UTR_info, UTR_seq in UTRs_dict.items():
        UTR_info_split=UTR_info.split('__')
        UTR_strand=UTR_info_split[2]
        UTR_name=UTR_info_split[3]
        
        if UTR_info in UTRs_BL21_dict:
            Reversal=UTRs_BL21_dict[UTR_info][2]
    
            if Reversal==0:
                if UTR_strand=='forward':
                    UTRs_BL21_dict[UTR_info].append('+')
                elif UTR_strand=='reverse':
                    UTRs_BL21_dict[UTR_info].append('-')
            elif Reversal==1:
                if UTR_strand=='forward':
                    UTRs_BL21_dict[UTR_info].append('-')
                elif UTR_strand=='reverse':
                    UTRs_BL21_dict[UTR_info].append('+')

    
    return UTRs_BL21_dict

UTRs_BL21_info_dict=parse_blast_results(Blast_results, UTRs_info_seq_dict)


#######
#Write BED file for RSeQC.
#######

def write_bed_and_broadPeak(UTRs_dict, deletions_data, pwd):
    print(len(UTRs_dict))
    
    UTRs_bed_file=open(pwd+"TUs_BL21.bed", 'w')
    UTRs_broadpeak_file=open(pwd+"TUs_BL21.broadPeak", 'w')
    UTRs_bed_file_del_cor=open(pwd+"TUs_BL21_del_cor.bed", 'w')
    UTRs_broadpeak_file_del_cor=open(pwd+"TUs_BL21_del_cor.broadPeak", 'w')    
    
    for UTR_name_info, UTR_info in UTRs_dict.items():
        UTR_name_info_split=UTR_name_info.split('__')
        Chrom_ID='NC_012971.2'
        UTR_chromStart=UTR_info[0]
        UTR_chromEnd=UTR_info[1]
        UTR_name=UTR_name_info_split[3] 
        UTR_score=100
        UTR_strand=UTR_info[3]
        UTR_thickStart=UTR_chromStart-1
        UTR_thickEnd=UTR_chromEnd
        UTR_itemRgb='255,0,0'
        UTR_blockCount=1
        UTR_blockSizes=UTR_chromEnd-UTR_chromStart+1
        UTR_blockStarts='0,'
        
        UTRs_bed_file.write(f'{Chrom_ID}\t{UTR_chromStart}\t{UTR_chromEnd}\t{UTR_name}\t{UTR_score}\t{UTR_strand}\t{UTR_thickStart}\t{UTR_thickEnd}\t{UTR_itemRgb}\t{UTR_blockCount}\t{UTR_blockSizes}\t{UTR_blockStarts}\n')
        
        UTR_signalValue=100
        UTR_pValue=-1
        UTR_qValue=-1
        
        UTRs_broadpeak_file.write(f'{Chrom_ID}\t{UTR_chromStart}\t{UTR_chromEnd}\t{UTR_name}\t{UTR_score}\t{UTR_strand}\t{UTR_signalValue}\t{UTR_pValue}\t{UTR_qValue}\n')
        
        Test_delition=0
        for deletion in deletions_data:
            if ((deletion[1]>=UTR_chromStart) & (UTR_chromStart>=deletion[0])) & ((deletion[1]>=UTR_chromEnd) & (UTR_chromEnd>=deletion[0])):
                Test_delition=1
        if Test_delition==0:
            UTRs_bed_file_del_cor.write(f'{Chrom_ID}\t{UTR_chromStart}\t{UTR_chromEnd}\t{UTR_name}\t{UTR_score}\t{UTR_strand}\t{UTR_thickStart}\t{UTR_thickEnd}\t{UTR_itemRgb}\t{UTR_blockCount}\t{UTR_blockSizes}\t{UTR_blockStarts}\n')
            UTRs_broadpeak_file_del_cor.write(f'{Chrom_ID}\t{UTR_chromStart}\t{UTR_chromEnd}\t{UTR_name}\t{UTR_score}\t{UTR_strand}\t{UTR_signalValue}\t{UTR_pValue}\t{UTR_qValue}\n')
        
    UTRs_bed_file.close()
    UTRs_broadpeak_file.close()
    UTRs_bed_file_del_cor.close()
    UTRs_broadpeak_file_del_cor.close()    
        
    return

write_bed_and_broadPeak(UTRs_BL21_info_dict, Deletions_data, PWD)