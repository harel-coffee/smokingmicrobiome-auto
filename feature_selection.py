#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:40:33 2020

@author: diego
"""

import os
import subprocess
import pandas as pd
import numpy as np
from collections import Counter
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelBinarizer
from scipy.stats import mannwhitneyu
from library.MLfunctions import get_random_samples


def fdr(p_vals):
    '''
    Computes Benjamini-Hochberg pvalue correction
    Parameters
    ----------
    p_vals : array
        list of p-values.

    Returns
    -------
    fdr : array
        list of corrected p-values.

    '''
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def fwdr(p_vals):
    '''
    Computes Bonferroni pvalue correction
    Parameters
    ----------
    p_vals : array
        list of p-values.

    Returns
    -------
    fdr : array
        list of corrected p-values.

    '''
    fwdr = p_vals * len(p_vals)
    fwdr[fwdr > 1] = 1

    return fwdr



''' 
def get_random_samples(df_metadata, target):
            
    t = 0.2
    phenotype = list(set(df_metadata[target]))
    random_samples = []
    np.random.seed()  
    for p in phenotype:        
        tmp = df_metadata.loc[df_metadata[target] == p]        
        for j in np.random.choice(tmp, int(t*len(tmp)), replace=False):
            random_samples.append(j)
    return random_samples
'''



def total_sum(df_pos):    
    for i, j in df_pos.iterrows():
        df_pos.loc[i] = j/sum(j)
    return df_pos
    

def create_tmp_dirs(folder):
    
    if not os.path.isdir(folder):
        cmd = 'mkdir '+folder
        subprocess.call(cmd, shell=True)        
        return True
    

if __name__ == "__main__":
        
    dirpath = os.sys.argv[1] # output folder
    dirpath = dirpath + "feature_selection/"
    filename_pos = os.sys.argv[2] # "data/species_intersect.csv"    
    filename_pheno = os.sys.argv[3] # "data/metadata.csv"
    target = os.sys.argv[4] # target
    id_samples = os.sys.argv[5] # target
    
    #dirpath = "data/results/"
    #dirpath = dirpath + "feature_selection/"
    #filename_pos = "data/species_intersect.csv"
    #filename_pheno = "data/metadata.csv"        
    create_tmp_dirs(dirpath)    
    n = 10  # number of iterations
    perc = 0.9 # percentage to allow taxas overall e.g. 0.9 = 90%        
    df_dict_taxa = pd.DataFrame()
    all_x_train = []
    all_y_train = []
    
    for i in range(n):        
        df_metadata = pd.read_csv(filename_pheno, index_col = 0)         
        
        #df_metadata.replace("former","never", inplace=True) #replace former
        
        df_pos = pd.read_csv(filename_pos, index_col=0)        
        random_samples_test = get_random_samples(df_metadata, target, id_samples)        
        df_pos = df_pos.loc[~df_pos.index.isin(random_samples_test)] # ignore test samples ~    
        samples = np.intersect1d(df_metadata.SRR, df_pos.index )
        df_metadata = df_metadata.loc[df_metadata[id_samples].isin(df_pos.index)]
        df_pos = df_pos.loc[df_metadata.SRR]                          
        lb = LabelBinarizer()
        #y_labels = df_metadata["smoking_status"].values
        y_labels = df_metadata[target].values
        y = lb.fit_transform(y_labels).flatten()                
        X = df_pos.values
        cv_splits           = 5 # outer loop
        fold                = 0           
        kfold               = StratifiedKFold(n_splits = cv_splits, shuffle = True)            
        for index_train, index_test in kfold.split(X, y):            
            fold += 1
            #print(f"Fold #{fold}")
            X_train = X[index_train]
            y_train = y[index_train]            
            X_test = X[index_test]
            y_test = y[index_test]                              
            X_train = total_sum(pd.DataFrame(X_train)).values                                
            X_test = total_sum(pd.DataFrame(X_test)).values                  
            df_x_train = pd.DataFrame(X_train)
            df_x_train.index = y_train
            df_x_train.columns = df_pos.columns
            all_x_train.append(df_x_train)
            all_y_train.append(y_train)
            dict_taxa_sig = {}
            for taxa in df_x_train.columns:
                group1 = df_x_train[taxa].loc[0]
                group2 = df_x_train[taxa].loc[1]
                u, p_value = mannwhitneyu(group1, group2)
                #print("two-sample wilcoxon-test", p_value)
                dict_taxa_sig[taxa] = p_value
            tmp = pd.DataFrame([dict_taxa_sig])  
            df_dict_taxa = pd.concat([df_dict_taxa, tmp])        
    df_all_x_train = pd.concat(all_x_train)
    df_dict_taxa_pvals = pd.DataFrame(list(map(fdr, df_dict_taxa.values))) #BH    
    df_dict_taxa_pvals.columns = df_pos.columns
    threshold = 0.05 # p-value threshold after correction
    dict_taxa_sig_corrected = {}
    for taxa in df_dict_taxa_pvals.columns:
        dict_taxa_sig_corrected[taxa] = sum(df_dict_taxa_pvals[taxa] < threshold)
    sigs = sorted(dict_taxa_sig_corrected.values())
    Counter(sigs)
    threshold_folds = int((fold*n)*perc)
    # consensus that appear at least threshold_folds from all splits
    sigs_taxas = {k:v for (k,v) in dict_taxa_sig_corrected.items() if v >= threshold_folds} 
    consensus_taxas = list(sigs_taxas.keys())
    #df_all_x_train = df_all_x_train[consensus_taxas]
    df_dict_taxa #uncorrected
    df_dict_taxa_pvals #corrected

    dict(sorted(sigs_taxas.items(), key=lambda item: item[1]))
    dict(sorted(dict_taxa_sig_corrected.items(), key=lambda item: item[1]))
    #print(df_all_x_train.shape)
    #df_all_x_train.to_csv(dirpath + "all_train_folds_feature_selection.csv")
    df_dict_taxa.to_csv(dirpath + "taxas_wilcoxon_uncorrected.csv")
    df_dict_taxa_pvals.to_csv(dirpath + "taxas_wilcoxon_BH_corrected.csv")

    #df_pos = pd.read_csv(filename_pos, index_col=0)
    #df_pos_consensus = df_pos[consensus_taxas]
    pd.DataFrame(consensus_taxas).to_csv(dirpath + "features.txt", index = None, header = False)
    #df_pos_consensus.to_csv(dirpath + "species_feature_selection_consensus.csv")
    