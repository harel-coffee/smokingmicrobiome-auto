#!/usr/bin/env python

# Copyright (C) 2021 Diego Montiel Gonzalez
# Erasmus Medical Center
# Department of Genetic Identification
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# TADA preprocessing

import os
import subprocess
import pandas as pd
from argparse   import ArgumentParser


def get_arguments():

    parser = ArgumentParser()    
    parser.add_argument("-c", "--consensus",  dest="consensus",
            help="Species_OTU_seqs_consensus.csv", metavar="FILE", required=True)    
    parser.add_argument("-s", "--species",  dest="species",
            help="species_intersect.csv", metavar="FILE", required=True)        
    parser.add_argument("-o", "--tada_out",  dest="tada_out",
            help="data/TADA/", metavar="DIR", required=True)                
    args = parser.parse_args()    
    return args


if __name__ == "__main__":

    args = get_arguments() 
    
    filename_ASV_consensus = args.consensus
    filename_species = args.species
    tada_output = args.tada_out

    #filename_ASV_consensus = os.sys.argv[1] # Species_OTU_seqs_consensus
    #filename_species = os.sys.argv[2] # S1_S5_species_intersect
    #tada_output = os.sys.argv[3] # TADA dirpath

    df_pos = pd.read_csv(filename_species, index_col=0)    
    df_consensus = pd.read_csv(filename_ASV_consensus, index_col=0)
    df_consensus = df_consensus.loc[df_pos.columns]
    df_consensus.index == df_pos.columns

    c = 1    
    file_fasta = tada_output + "/taxa_species.fasta"
    file = open(file_fasta,"w")
    list_cols = []

    for i in df_consensus["sequence"].values:        
        
        #file.write(">" + df_consensus.index[c])
        file.write(">"+ str(c))
        file.write("\n")
        file.write(str(i))
        file.write("\n")
        c+=1

    file.close()      
    file_mafft = file_fasta.split(".")[:-1][0] + "_mafft.fasta"

    if not os.path.exists(file_mafft):
        #### MAFFT    
        cmd = "mafft --thread -1 --auto {} > {}".format(file_fasta, file_mafft)
        subprocess.check_output(cmd, shell=True)    
        #### FastTree
        file_tree = file_fasta.split(".")[:-1][0] + "_tree.nwk"
        cmd = "FastTree -gtr -nt {} > {}".format(file_mafft, file_tree)
        subprocess.check_output(cmd, shell=True)

    print(os.listdir(tada_output))
