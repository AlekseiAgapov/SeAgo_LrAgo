###############################################
##Dmitry Sutormin, 2022##
##Ago-Seq analysis##

####
#Convert coverage depth tracks to RPKM tracks.
####

###############################################


#Path to the input coverage files.
PWD='C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Ago_Seq\\'
filein_path_dict={'1'  : PWD + "Syn_Lro_libraries_data_wig_polyshed\Se_5.5_nocip_rep_1_minus.wig",
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

#Path to the output RPKM files.
fileout_path_dict={'1'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_nocip_rep_1_minus_FPKM.wig",
                   '2'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_nocip_rep_1_plus_FPKM.wig",
                   '3'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_nocip_rep_2_minus_FPKM.wig",
                   '4'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_nocip_rep_2_plus_FPKM.wig",
                   '5'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_cip_rep_1_minus_FPKM.wig",
                   '6'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_cip_rep_1_plus_FPKM.wig",
                   '7'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_cip_rep_2_minus_FPKM.wig",
                   '8'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_5.5_cip_rep_2_plus_FPKM.wig",
                   '9'  : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_nocip_rep_1_minus_FPKM.wig",
                   '10' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_nocip_rep_1_plus_FPKM.wig",
                   '11' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_nocip_rep_2_minus_FPKM.wig",
                   '12' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_nocip_rep_2_plus_FPKM.wig",
                   '13' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_cip_rep_1_minus_FPKM.wig",
                   '14' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_cip_rep_1_plus_FPKM.wig",
                   '15' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_cip_rep_2_minus_FPKM.wig",
                   '16' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Se_12.5_cip_rep_2_plus_FPKM.wig",
                   '17' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_nocip_rep_1_minus_FPKM.wig",
                   '18' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_nocip_rep_1_plus_FPKM.wig",
                   '19' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_nocip_rep_2_minus_FPKM.wig",
                   '20' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_nocip_rep_2_plus_FPKM.wig",
                   '21' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_cip_rep_1_minus_FPKM.wig",
                   '22' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_cip_rep_1_plus_FPKM.wig",
                   '23' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_cip_rep_2_minus_FPKM.wig",
                   '24' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_5.5_cip_rep_2_plus_FPKM.wig",
                   '25' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_nocip_rep_1_minus_FPKM.wig",
                   '26' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_nocip_rep_1_plus_FPKM.wig",
                   '27' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_nocip_rep_2_minus_FPKM.wig",
                   '28' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_nocip_rep_2_plus_FPKM.wig",
                   '29' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_cip_rep_1_minus_FPKM.wig",
                   '30' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_cip_rep_1_plus_FPKM.wig",
                   '31' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_cip_rep_2_minus_FPKM.wig",
                   '32' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\Lr_12.5_cip_rep_2_plus_FPKM.wig",
                    }


#Path to the reads numbers data.
read_num_dict={'1'  : 368873,
               '2'  : 368873,
               '3'  : 488441,
               '4'  : 488441,
               '5'  : 1509243,
               '6'  : 1509243,
               '7'  : 1038304,
               '8'  : 1038304,
               '9'  : 791799,
               '10' : 791799,
               '11' : 675406,
               '12' : 675406,
               '13' : 2152542,
               '14' : 2152542,
               '15' : 2920406,
               '16' : 2920406,
               '17' : 1266921,
               '18' : 1266921,
               '19' : 1984106,
               '20' : 1984106,
               '21' : 1290785,
               '22' : 1290785,
               '23' : 1265045,
               '24' : 1265045,
               '25' : 1196750,
               '26' : 1196750,
               '27' : 1003480,
               '28' : 1003480,
               '29' : 2974923,
               '30' : 2974923,
               '31' : 1556966,
               '32' : 1556966,
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


def read_and_convert_to_RPKM(filein_path_dict, fileout_path_dict, read_num_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        scaling_factor=float(1)*1000000*1000/int(read_num_dict[sample_name])
        
        for line in filein:
            line=line.rstrip()
            if line[0] in ['t']:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif line[0] in ['f']:
                continue
            else:
                fileout.write(str(int(line)*scaling_factor)+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert_to_RPKM(filein_path_dict, fileout_path_dict, read_num_dict, name_dict, Chromosome_name, Auto_or_manual)


