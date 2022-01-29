#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Copyright (C) 2021 Diego Montiel Gonzalez
# Erasmus Medical Center
# Department of Genetic Identification
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# NestedCV pipeline for Smoking Microbiome sequencing data



import os
import subprocess
import warnings
import library.MLfunctions as func

import pandas as pd
import numpy as np

from collections import Counter
from argparse   import ArgumentParser


from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelBinarizer
from sklearn.metrics import make_scorer, classification_report

import matplotlib.pyplot as plt
import seaborn as sns
import joblib


def get_arguments():

    parser = ArgumentParser()    

    parser.add_argument("-i", "--input", dest = "filename_pos", required = True,  
        help="table taxa species e.g. data/species_intersect.csv. Taxa as columns, sample names as rows", metavar="FILE")

    parser.add_argument("-o", "--output", dest = "output_dir", required = True,
        help = "output folder name e.g. data/output/", metavar = "PATH")            

    parser.add_argument("-m", "--metadata", dest = "filename_pheno", required = True,
        help="table with metadata in Tab-separated values TSV e.g. data/metadata.tsv", metavar="FILE")

    parser.add_argument("-target", "--target", dest = "target", required = True,
        help="Target phenotype name based on the column of metadata e.g. smoking", metavar="STRING")

    parser.add_argument("-id", "--id", dest = "id", required = True,
        help="Column id sample names on the column of metadata e.g. SRR or Samples", metavar="STRING")

    parser.add_argument("-t", "--threads",  dest="threads", default = 4, type = int,
        help="threads to be used during Nested Cross-Validation e.g. 4", 
        metavar="INT", required=False) 
    
    parser.add_argument("-r", "--root-tree", dest="tree_newick", required = False, default = False,
        help="Rooted tree in newick format from ASV consensus e.g. data/TADA/taxa_species_tree.nwk", 
        metavar="FILE")
    
    parser.add_argument("-ml", "--ml-model", required = False, dest = "model", default = False,
        help="Machine learning models, if more than one is selected implements NestedCV and\
        reports output and prediction estimations using MCC and AUC,\
        otherwise you can select one to export as final model (also needs to be selected an augmented data type)\
        [LR, KNN, DT, SVMR, SVML, RF, XG] e.g. LR,KNN,RF or e.g. SVML (will generate one model and export)")
    
    parser.add_argument("-a", "--augmented", required = False, dest = "augmented", default = False,
        help = "Data type augmented technique, if no type is selected uses all data types\
        [DEFAULT, ADASYN_over, ADASYN_both, SMOTE_both,  SMOTE_over, TADA] ")
    
    parser.add_argument("-iter", "--iterations", dest = "iterations", required = False,
        help="number of iterations [DEFAULT=10]", metavar="INT", type = int, default = 10)

            
    args = parser.parse_args()
    return args


def create_tmp_dirs(folder):
    
    if not os.path.isdir(folder):
        cmd = 'mkdir '+folder
        subprocess.call(cmd, shell=True)        
        return True
    

if __name__ == "__main__":
    
    args = get_arguments()
    filename_pos = args.filename_pos
    output_dir = args.output_dir + "/"
    filename_pheno = args.filename_pheno
    
    ## optionals    
    tree_newick = args.tree_newick    
    methods = args.augmented
    estimators = args.model
    target = args.target
    id_samples = args.id
    
    threads = args.threads
    iterations = args.iterations
    
    if methods == False:
        methods = ["DEFAULT", "ADASYN_over", "ADASYN_both", "SMOTE_both",  "SMOTE_over"]
        if tree_newick:
            methods.append("TADA")
    else:
        methods = methods.split(",")        
    if estimators == False:
        estimators = ["LR", "KNN", "DT", "SVMR", "SVML", "RF", "XG"]
    else:
        estimators = estimators.split(",")
        
    print("-- Running --")
    
    print("Input taxa species: {}".format(filename_pos))
    print("Output directory: {}".format(output_dir))
    print("Metadata: {}".format(filename_pheno))    
    print("You selected the following model(s):")    
    print(" ".join(estimators))
    print("You selected the following augmentation method(s):")    
    print(" ".join(methods))    
    print("Number of threads: {}".format(threads))
    
    flag_export_model = False
    # If true export models and test set from nestedCV
    if len(estimators) == 1 == len(methods):        
        flag_export_model = True                   
    create_tmp_dirs(output_dir)        
    mcc_test_list = []
    auc_test_list = []    
    pd_all_metrics = pd.DataFrame()    
    for i in range(iterations):
        print(i+1)
        df_metadata = pd.read_csv(filename_pheno, index_col = 0, sep="\t")
        df_pos = pd.read_csv(filename_pos, index_col=0)    
        taxas = df_pos.columns
        df_pos.columns = list(range(1, len(df_pos.columns)+1))
        print(df_metadata.head())
        random_samples_test = func.get_random_samples(df_metadata, target, id_samples)        
        df_pos = df_pos.loc[~df_pos.index.isin(random_samples_test)] # ignore test samples ~            
        #print(random_samples_test)
        df_metadata = df_metadata.loc[df_metadata[id_samples].isin(df_pos.index)]
        df_pos = df_pos.loc[df_metadata[id_samples]]
        lb = LabelBinarizer()
        y_labels = df_metadata[target].values    
            
        y = lb.fit_transform(y_labels).flatten()                
        X = df_pos.values
            
        for choose in methods:
            
            for estimator in estimators:        
                pd_model            = pd.DataFrame()
                list_labels         = []                
                auc_kfold           = []
                mcc_kfold           = []           
                val_loss_list       = []
                val_acc_list        = []    
                y_probs_max_list    = []
                accuracy_list       = []
                cv_splits           = 5 # outer loop            
                inner_cv_splits     = 2            
                fold                = 0 
                y_probs             = np.array(0)   
                list_names_kfold    = []   
                df_activations      = pd.DataFrame()                
                kfold               = StratifiedKFold(n_splits=cv_splits, shuffle=True)
                
                print("############################################")                
                print("Augmentation: {} ML Model: {}".format(choose, estimator))
                print("############################################")
                
                metrics_models = {}
                scores_mcc = []
                scores_auc = []
                
                for index_train, index_test in kfold.split(X, y):          
                    print("------------------")
                    fold+=1
                    print(f"Fold #{fold}")
                    X_train = X[index_train]
                    y_train = y[index_train]            
                    X_test = X[index_test]
                    y_test = y[index_test]                    
                    if choose == "DEFAULT":            
                        pass
                    if choose == "TADA":                    
                        obs_ids = [str(i) for i in df_pos.columns]                       
                        X_train, y_train = func.tada_augmentation(X_train, y_train,  obs_ids, fold, tree_newick)                
                    elif  choose == "ADASYN_over":                
                        ada = ADASYN(sampling_strategy='minority', n_neighbors = 5, n_jobs = threads)
                        X_train, y_train = ada.fit_resample(X_train,y_train)                
                    elif choose == "SMOTE_both":                        
                        X_train, y_train = func.smote_both(X_train, y_train)     
                    
                    elif choose == "ADASYN_both":                        
                        X_train, y_train = func.adasyn_both(X_train, y_train)     
                    
                    elif choose == "SMOTE_over":                      
                        X_train, y_train = SMOTE().fit_resample(X_train, y_train)                                                        
                    X_train = func.total_sum(pd.DataFrame(X_train)).values                        
                    X_test = func.total_sum(pd.DataFrame(X_test)).values                                                        
                    pd_model = pd.concat([pd_model, pd.DataFrame(X_train)])    
                    
                    list_labels.append(y_train)        
                    cv = StratifiedKFold(n_splits = inner_cv_splits, shuffle=True)                    
                    if "LR" in estimator:
                        clf = func.LogisticRegression()     
                        param_grid = func.GridParams().getLogistic()                    

                    elif "KNN" in estimator:            
                        clf = func.KNeighborsClassifier()                      
                        param_grid = func.GridParams().getKNN()
                                
                    elif "SVMR" in estimator:
                        clf = func.SVC(probability = True)
                        param_grid = func.GridParams().getSVMRadial()
                        
                    elif "SVML" in estimator:            
                        clf = func.SVC(probability = True)
                        param_grid = func.GridParams().getSVMLinear()        
                    
                    elif "DT" in estimator:
                        clf = func.DecisionTreeClassifier()
                        param_grid = func.GridParams().getDT()                                       
                                        
                    elif "RF" in estimator:
                        clf = func.RandomForestClassifier()
                        param_grid = func.GridParams().getRF()
                        
                    elif "XG" in estimator:                   
                        clf = func.XGBClassifier(objective = 'binary:logistic')                
                        param_grid = func.GridParams().getXG()                                            
                    best_algo = func.RandomizedSearchCV(estimator = clf, 
                                        param_distributions = param_grid, cv = cv, 
                                        scoring = make_scorer(func.calculate_mcc),
                                        verbose = 1, n_jobs = threads, n_iter=1000) 
                    best_algo.fit(X_train,y_train)
                                                    
                    if flag_export_model:
                        model_path = output_dir + "/Model/"
                        create_tmp_dirs(model_path)
                        joblib.dump(best_algo, model_path + "model_" + str(i) + "_"
                                    + str(fold) + ".model")                                                                
                    y_pred = best_algo.predict(X_test)
                    y_probs = best_algo.predict_proba(X_test)        
                    test_mcc, test_auc, test_cnf_matrix = func.get_metrics_classification(y_test, y_pred, y_probs)    
                    test_auc = test_auc[0]                                                
                    out_report = classification_report(y_test,y_pred, output_dict = True,
                                                    target_names = ["smoker", "non-smoker"])                        
                    if estimator + "_mcc" in metrics_models:
                        metrics_models[estimator + "_mcc"].append(test_mcc)
                    else:
                        metrics_models[estimator + "_mcc"] = [test_mcc]                            
                    if estimator + "_auc" in metrics_models:
                        metrics_models[estimator + "_auc"].append(test_auc)
                    else:
                        metrics_models[estimator + "_auc"] = [test_auc]         
                    if estimator + "_f1" in metrics_models:
                        metrics_models[estimator + "_f1"].append(out_report["macro avg"]["f1-score"])
                    else:
                        metrics_models[estimator + "_f1"] = [out_report["macro avg"]["f1-score"]]                     
                                    
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')    
                    pd_all_metrics = pd.concat([pd_all_metrics, func.process_metrics_df(metrics_models, choose)])        
        if flag_export_model: 
            mcc_test, auc_test = func.get_test_set(filename_pheno, filename_pos, random_samples_test, model_path, i, id_samples, target)
            mcc_test_list.append(mcc_test)
            auc_test_list.append(auc_test[0])    
    
    #<-------------------------------------------------------------
    result = " Mean NestedCV train"
    result += "\n"
    for type_ in set(pd_all_metrics["type"]):            
        result += "Augmentation technique: {}".format(type_)
        #result += type_
        result += "\n"                
        for model in  set(pd_all_metrics["model"]):
            result +=  "Model: {}".format(model)
            #result +=  model
            result += "\n"                
            mcc_mean = pd_all_metrics.loc[(pd_all_metrics["model"] == model) &
                                (pd_all_metrics["type"] == type_)]["mcc"].mean()
            mcc_std = pd_all_metrics.loc[(pd_all_metrics["model"] == model) &
                                (pd_all_metrics["type"] == type_)]["mcc"].std()                    
            auc_mean = pd_all_metrics.loc[(pd_all_metrics["model"] == model) &
                                (pd_all_metrics["type"] == type_)]["auc"].mean()
            auc_std = pd_all_metrics.loc[(pd_all_metrics["model"] == model) &
                                (pd_all_metrics["type"] == type_)]["auc"].std()
            result += "MCC mean: {} Std: {} +/-".format(round(mcc_mean,3), round(mcc_std,3))
            result += "\n"                
            result += "AUC mean: {} Std: {} +/-".format(round(auc_mean,3), round(auc_std,3))
            result += "\n"
            result += "-" * 50
            result += "\n"
        result += "-" * 50
        result += "\n"    
    
    if not flag_export_model:
        file_mean_results = "mean_results_train.txt"
        pd_all_metrics.to_csv(output_dir + "NestedCV_train_results.csv", sep="\t")        
        plt.figure(figsize=(20, 15))
        x = sns.boxplot(x="mcc", y="model", hue="type", data=pd_all_metrics, 
                    dodge=True, linewidth=1.2, palette = "Set3")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(output_dir + "mcc_boxplot.png")
        
        plt.figure(figsize=(20, 15))
        x = sns.boxplot(x="auc", y="model", hue="type", data=pd_all_metrics, 
                    dodge=True, linewidth=1.2, palette = "Set3")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(output_dir + "auc_boxplot.png")
    
    if flag_export_model:
        file_mean_results = "mean_results_test.txt"
        pd_all_metrics.to_csv(output_dir + "NestedCV_test_results.csv", sep="\t")        
        result += "\nTest set"
        result += "\n"
        result += 'MCC: {} std +/- {}'.format(np.mean(mcc_test_list), np.std(mcc_test_list))
        result += "\n"
        result += 'AUC: {} std +/- {}'.format(np.mean(auc_test_list), np.std(auc_test_list))
        result += "\n"        
        
    with open(output_dir + file_mean_results, 'w+') as f:
        f.write(result)
   
    print("-- NestedCV pipeline finished---")
    print("Check {} for output generated!".format(output_dir))

