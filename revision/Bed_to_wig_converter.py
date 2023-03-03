###############################################
##Dmitry Sutormin, 2023##
##Ago-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Ago_Seq\\'
filein_path_dict={'1'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\01_minus_cov.tsv",
                  '2'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\01_plus_cov.tsv",
                  '3'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\02_minus_cov.tsv",
                  '4'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\02_plus_cov.tsv",
                  '5'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\03_minus_cov.tsv",
                  '6'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\03_plus_cov.tsv",
                  '7'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\04_minus_cov.tsv",
                  '8'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\04_plus_cov.tsv",
                  '9'  : PWD + "Syn_Lro_libraries_data_bed_polyshed\\05_minus_cov.tsv",
                  '10' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\05_plus_cov.tsv",
                  '11' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\06_minus_cov.tsv",
                  '12' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\06_plus_cov.tsv",
                  '13' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\07_minus_cov.tsv",
                  '14' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\07_plus_cov.tsv",   
                  '15' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\08_minus_cov.tsv",
                  '16' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\08_plus_cov.tsv",
                  '17' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\09_minus_cov.tsv",
                  '18' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\09_plus_cov.tsv",
                  '19' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\10_minus_cov.tsv",
                  '20' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\10_plus_cov.tsv",
                  '21' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\11_minus_cov.tsv",
                  '22' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\11_plus_cov.tsv",
                  '23' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\12_minus_cov.tsv",
                  '24' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\12_plus_cov.tsv",
                  '25' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\13_minus_cov.tsv",
                  '26' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\13_plus_cov.tsv",
                  '27' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\14_minus_cov.tsv",
                  '28' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\14_plus_cov.tsv", 
                  '29' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\15_minus_cov.tsv",
                  '30' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\15_plus_cov.tsv",
                  '31' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\16_minus_cov.tsv",
                  '32' : PWD + "Syn_Lro_libraries_data_bed_polyshed\\16_plus_cov.tsv",                    
                  }

#Path to the output file.
fileout_path_dict={'1'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_nocip_rep_1_minus.wig",
                   '2'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_nocip_rep_1_plus.wig",
                   '3'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_nocip_rep_2_minus.wig",
                   '4'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_nocip_rep_2_plus.wig",
                   '5'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_cip_rep_1_minus.wig",
                   '6'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_cip_rep_1_plus.wig",
                   '7'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_cip_rep_2_minus.wig",
                   '8'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_cip_rep_2_plus.wig",
                   '9'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_nocip_rep_1_minus.wig",
                   '10' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_nocip_rep_1_plus.wig",
                   '11' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_nocip_rep_2_minus.wig",
                   '12' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_nocip_rep_2_plus.wig",
                   '13' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_cip_rep_1_minus.wig",
                   '14' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_cip_rep_1_plus.wig",
                   '15' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_cip_rep_2_minus.wig",
                   '16' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_12.5_cip_rep_2_plus.wig",
                   '17' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_nocip_rep_1_minus.wig",
                   '18' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_nocip_rep_1_plus.wig",
                   '19' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_nocip_rep_2_minus.wig",
                   '20' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_nocip_rep_2_plus.wig",
                   '21' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_cip_rep_1_minus.wig",
                   '22' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_cip_rep_1_plus.wig",
                   '23' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_cip_rep_2_minus.wig",
                   '24' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_5.5_cip_rep_2_plus.wig",
                   '25' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_nocip_rep_1_minus.wig",
                   '26' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_nocip_rep_1_plus.wig",
                   '27' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_nocip_rep_2_minus.wig",
                   '28' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_nocip_rep_2_plus.wig",
                   '29' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_cip_rep_1_minus.wig",
                   '30' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_cip_rep_1_plus.wig",
                   '31' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_cip_rep_2_minus.wig",
                   '32' : PWD + "Syn_Lro_libraries_data_wig_polyshed\Lr_12.5_cip_rep_2_plus.wig",
                    }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1'  : "Se_5.5_nocip_rep_1_minus.wig",
           '2'  : "Se_5.5_nocip_rep_1_plus.wig",
           '3'  : "Se_5.5_nocip_rep_2_minus.wig",
           '4'  : "Se_5.5_nocip_rep_2_plus.wig",
           '5'  : "Se_5.5_cip_rep_1_minus.wig",
           '6'  : "Se_5.5_cip_rep_1_plus.wig",
           '7'  : "Se_5.5_cip_rep_2_minus.wig",
           '8'  : "Se_5.5_cip_rep_2_plus.wig",
           '9'  : "Se_12.5_nocip_rep_1_minus.wig",
           '10' : "Se_12.5_nocip_rep_1_plus.wig",
           '11' : "Se_12.5_nocip_rep_2_minus.wig",
           '12' : "Se_12.5_nocip_rep_2_plus.wig",
           '13' : "Se_12.5_cip_rep_1_minus.wig",
           '14' : "Se_12.5_cip_rep_1_plus.wig",
           '15' : "Se_12.5_cip_rep_2_minus.wig",
           '16' : "Se_12.5_cip_rep_2_plus.wig",
           '17' : "Lr_5.5_nocip_rep_1_minus.wig",
           '18' : "Lr_5.5_nocip_rep_1_plus.wig",
           '19' : "Lr_5.5_nocip_rep_2_minus.wig",
           '20' : "Lr_5.5_nocip_rep_2_plus.wig",
           '21' : "Lr_5.5_cip_rep_1_minus.wig",
           '22' : "Lr_5.5_cip_rep_1_plus.wig",
           '23' : "Lr_5.5_cip_rep_2_minus.wig",
           '24' : "Lr_5.5_cip_rep_2_plus.wig",
           '25' : "Lr_12.5_nocip_rep_1_minus.wig",
           '26' : "Lr_12.5_nocip_rep_1_plus.wig",
           '27' : "Lr_12.5_nocip_rep_2_minus.wig",
           '28' : "Lr_12.5_nocip_rep_2_plus.wig",
           '29' : "Lr_12.5_cip_rep_1_minus.wig",
           '30' : "Lr_12.5_cip_rep_1_plus.wig",
           '31' : "Lr_12.5_cip_rep_2_minus.wig",
           '32' : "Lr_12.5_cip_rep_2_plus.wig",
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name='NC_012971.2'
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(1)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)