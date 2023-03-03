###############################################
##Dmitry Sutormin, 2022##
##ChIP-Seq analysis##

####
#The only purpose - to compute by-position average of a set of wig files.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm


#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Other\E_coli_pAgo_and_topo_activity\Ago_Seq\\"

#Dictionary of replicas 
#'Replica name' : 'Path to wig file'
Dict_of_replicas1={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_cip_rep_1_minus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_cip_rep_2_minus_FPKM.wig",
                   }

Dict_of_replicas2={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_cip_rep_1_plus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_cip_rep_2_plus_FPKM.wig",
                   }

Dict_of_replicas3={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_nocip_rep_1_minus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_nocip_rep_2_minus_FPKM.wig",
                   }

Dict_of_replicas4={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_nocip_rep_1_plus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_5.5_nocip_rep_2_plus_FPKM.wig",
                   }

Dict_of_replicas5={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_cip_rep_1_minus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_cip_rep_2_minus_FPKM.wig",
                   }

Dict_of_replicas6={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_cip_rep_1_plus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_cip_rep_2_plus_FPKM.wig",
                   }

Dict_of_replicas7={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_nocip_rep_1_minus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_nocip_rep_2_minus_FPKM.wig",
                   }

Dict_of_replicas8={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_nocip_rep_1_plus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Se_12.5_nocip_rep_2_plus_FPKM.wig",
                   }

Dict_of_replicas9={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_cip_rep_1_minus_FPKM.wig",
                   'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_cip_rep_2_minus_FPKM.wig",
                   }

Dict_of_replicas10={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_cip_rep_1_plus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_cip_rep_2_plus_FPKM.wig",
                    }

Dict_of_replicas11={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_nocip_rep_1_minus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_nocip_rep_2_minus_FPKM.wig",
                    }

Dict_of_replicas12={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_nocip_rep_1_plus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_5.5_nocip_rep_2_plus_FPKM.wig",
                    }

Dict_of_replicas13={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_cip_rep_1_minus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_cip_rep_2_minus_FPKM.wig",
                    }

Dict_of_replicas14={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_cip_rep_1_plus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_cip_rep_2_plus_FPKM.wig",
                    }

Dict_of_replicas15={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_nocip_rep_1_minus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_nocip_rep_2_minus_FPKM.wig",
                    }

Dict_of_replicas16={'Replic 1' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_nocip_rep_1_plus_FPKM.wig",
                    'Replic 2' : PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed\\Lr_12.5_nocip_rep_2_plus_FPKM.wig",
                    }

#ID or short description of the track (will be the name of a track in IGV).
name1='Se_5.5_cip_minus_FPKM_av'
name2='Se_5.5_cip_plus_FPKM_av'
name3='Se_5.5_nocip_minus_FPKM_av'
name4='Se_5.5_nocip_plus_FPKM_av'
name5='Se_12.5_cip_minus_FPKM_av'
name6='Se_12.5_cip_plus_FPKM_av'
name7='Se_12.5_nocip_minus_FPKM_av'
name8='Se_12.5_nocip_plus_FPKM_av'
name9='Lr_5.5_cip_minus_FPKM_av'
name10='Lr_5.5_cip_plus_FPKM_av'
name11='Lr_5.5_nocip_minus_FPKM_av'
name12='Lr_5.5_nocip_plus_FPKM_av'
name13='Lr_12.5_cip_minus_FPKM_av'
name14='Lr_12.5_cip_plus_FPKM_av'
name15='Lr_12.5_nocip_minus_FPKM_av'
name16='Lr_12.5_nocip_plus_FPKM_av'
#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name='NC_012971.2'
#Output path for the final file.
PWD_out=PWD + "Syn_Lro_libraries_data_wig_RPKM_polyshed_av\\"



#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return

#########
##Write file with avaraged data.
#########

def write_wig_file(dict_of_replicas, average_file_path, name, Chromosome_name):
    
    average_out=open(average_file_path, 'w')
    average_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    
    for i in range(len(dict_of_replicas[list(dict_of_replicas.keys())[0]])):
        av_data_position=[]
        for replica_name, replica_data in dict_of_replicas.items():
            av_data_position.append(replica_data[i])
        average_out.write(str(np.mean(av_data_position))+'\n')
    
    average_out.close()
    return

#########
##Wrapper function.
#########

def Wrapper_wig_av(dict_of_replicas_path, name, pwd_out, chromosome_name):
    
    #Contains data of all replicas in separate arrays.
    dict_of_replicas={}
    for replica_name, replica_path in dict_of_replicas_path.items():
        dict_of_replicas[replica_name]=wig_parsing(replica_path)   
        
    #Calculate correlation between replicates and draw heatmap.  
    outpath=pwd_out + name + "_correlation_matrix.png"
    correlation_matrix(pd.DataFrame(dict_of_replicas), 'pearson', 'Correlation of biological replicas', outpath)
    
    #Average tracks and write resultant .wig file.
    average_file_path=pwd_out + name + ".wig"
    write_wig_file(dict_of_replicas, average_file_path, name, chromosome_name)
    
    return

Wrapper_wig_av(Dict_of_replicas1, name1, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas2, name2, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas3, name3, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas4, name4, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas5, name5, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas6, name6, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas7, name7, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas8, name8, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas9, name9, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas10, name10, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas11, name11, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas12, name12, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas13, name13, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas14, name14, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas15, name15, PWD_out, Chromosome_name)
Wrapper_wig_av(Dict_of_replicas16, name16, PWD_out, Chromosome_name)