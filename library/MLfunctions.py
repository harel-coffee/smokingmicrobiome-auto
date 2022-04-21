#!/usr/bin/env python

import os
import subprocess
import warnings
from collections import Counter
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import OneHotEncoder


from sklearn.metrics import matthews_corrcoef, make_scorer, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_fscore_support
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.tree import DecisionTreeClassifier

from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.under_sampling import RandomUnderSampler

from imblearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold


import matplotlib.pyplot as plt
from scipy import interp

from biom.table import Table
from biom import load_table
import biom
import joblib


class GridParams:
    
    # Randomized
    gamma_range = np.linspace(0.0001, 10, 100)
    c_range = np.linspace(1, 10, 1000)
        
    def getKNN(self):
        return {'n_neighbors': list(range(5,101)), 'weights': ["uniform", "distance"], 
                'algorithm': ["auto", "kd_tree", "ball_tree", "brute"], 
                'leaf_size': list(range(30,61)), 'p':[1,2]}    
    
    def getLogistic(self):
        return {'penalty' : ['l1'], 'C' : np.logspace(-4, 4, 20), 
                        'max_iter' : [1000], 'solver' : ['liblinear']}                                
    
    def getSVMRadial(self):
        return {'C': self.c_range, 'gamma': self.gamma_range, 'kernel': ["rbf"]}  
    
    def getSVMLinear(self):
        return {'C': self.c_range, 'kernel': ["linear"]}  
            
    def getDT(self):
        return {"criterion": ["gini","entropy"],
                "max_depth": range(1,21),
                "min_samples_split": range(2,21),
                "min_samples_leaf": range(1,10),
                "max_features": ["auto", "sqrt", "log2"]
                }
    
    def getRF(self):                
        return {    
            'n_estimators' : range(100, 1001, 50),
            'max_features': ["auto", "sqrt"],
            'max_depth' : range(10,101, 10),    
            'min_samples_split': range(2,11),
            'min_samples_leaf': range(1,5),
            'bootstrap': [True, False]
            }
     
    def getXG(self):         
        return {
                'n_estimators': range(100, 1001, 50),
                'learning_rate': np.linspace(0.05, 3, 10),
                'reg_alpha': np.linspace(0, 1, 10),
                'reg_lambda': np.linspace(0, 1, 10),
                'min_child_weight': range(1, 11),
                'max_depth': range(3, 11),
                'gamma' : np.linspace(0, 0.2, 10),                
                'subsample': np.linspace(0.5, 1, 10), 
                'colsample_bytree': np.linspace(0.5, 1, 10),
                'boosting_type': ['gbdt']                
            }

def get_metrics_classification(y_test, y_pred, y_probs):
    """
    Returns: MCC score, AUCs, Confusion matrix
    
    Parameters
    ----------
    y_test : numpy array
        array of assign classes e.g [0,1,0,2...]
    y_pred : numpy array
        
        array with predicted label classes e.g. [0,1,0,2...]
    """
    
    mcc_score = calculate_mcc(y_test, y_pred)
    auc_score = calculate_roc_auc(y_test,y_probs)                
    cnf_matrix = confusion_matrix(y_test, y_pred)
    return mcc_score, auc_score, cnf_matrix
        

def get_random_samples(df_metadata, target,id_samples):
    
    t = 0.2
    phenotype = list(set(df_metadata[target]))
    random_samples = []
    np.random.seed()  
    for p in phenotype:                        
        tmp = df_metadata.loc[df_metadata[target] == p]        
        tmp_target = tmp[id_samples].values        
        for j in np.random.choice(tmp_target, int(t*len(tmp_target)), replace=False):
            random_samples.append(j)        
    return random_samples


def get_test_set(filename_pheno, filename_pos, random_samples_test, 
                                 model_path, i, id_samples, target):    
    
    df_test_pheno = pd.read_csv(filename_pheno, index_col = 0, sep="\t")            
    df_test_pos = pd.read_csv(filename_pos, index_col=0)    
    df_test_pos = df_test_pos.loc[df_test_pos.index.isin(random_samples_test)]    
    df_test_pheno = df_test_pheno.loc[df_test_pheno[id_samples].isin(df_test_pos.index)]
        
    df_test_pos = df_test_pos.loc[df_test_pheno[id_samples]]        
    df_test_pos = total_sum(df_test_pos)        
    y_test_labels = df_test_pheno[target].values            
    y_test = y_test_labels
    X_test = df_test_pos.values                       
    y_test = LabelBinarizer().fit_transform(y_test_labels)    
    if len(set(y_test_labels)) == 3:
        y_test = np.argmax(y_test, axis=1)
    elif len(set(y_test_labels)) == 2:
        y_test = y_test.flatten()        
    #random_samples_test = df_test_pos.index + "_"+y_test_labels + "_"+df_test_pheno["study"]
    random_samples_test = df_test_pos.index + "_" +y_test_labels   
    dirpath        = model_path
    models_dirpath = model_path
    outpath        = model_path
    y_probs_final = np.array(0)
    models_dirpath = os.walk(models_dirpath)
    count = 0
    for dirpath, dirnames, filenames in models_dirpath:
        for filename in [f for f in filenames if f.endswith(".model")]:
            path_model = os.path.join(dirpath, filename)            
            count+=1            
            model = joblib.load(path_model)            
            y_pred = model.predict(X_test)            
            y_probs_final = y_probs_final + model.predict_proba(X_test)
    y_probs_final = y_probs_final/count
    y_pred = np.argmax(y_probs_final, axis=1)    
    df_probs_final = pd.DataFrame(y_probs_final, index = random_samples_test, columns=["current","never"])    
    df_probs_final.to_csv(outpath + str(i) + "probs_test.txt")
    y_real = []
    for i in df_probs_final.index:
        if "current" in i:
            y_real.append(0)
        elif "former" in i:
            y_real.append(1)
        elif "never" in i:
            y_real.append(1)
    y_real = np.array(y_real)    
    mcc_score, auc_score, cnf_matrix = get_metrics_classification(y_real, y_pred, y_probs_final)        
    return mcc_score, auc_score
    

def calculate_mcc(y_test, y_pred):
        if len(set(y_pred)) > 1:            
            return matthews_corrcoef(y_test, y_pred)
        else:
            return 0
        

def total_sum(df_pos):    
    for i, j in df_pos.iterrows():
        df_pos.loc[i] = j/sum(j)
    return df_pos


def set_one_hot_encoded(target):

    onehot_encoder = OneHotEncoder(sparse=False, categories='auto')    
    integer_encoded = target.reshape(len(target), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)    
    
    return onehot_encoded
    

def calculate_roc_auc(y, y_score):
    """
    Parameters
    ----------
    y : numpy array
        array of assign classes e.g np[0,1,0,2...].
    y_score : numpy array
        matrix with assign probabilities.
    Returns
    Dict of RoC AUC values per class
    -------
    None.
    """
    onehot_encoded = set_one_hot_encoded(y)    
    n_classes = onehot_encoded.shape[1]            
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    
    for i in range(n_classes):
        fpr[i], tpr[i], thresholds = roc_curve(onehot_encoded[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])                
    
    return roc_auc


def create_tmp_dirs(folder):

    flag = True
    if os.path.isdir(folder):    
        while(flag):
            print("WARNING! File "+folder+" already exists, \nWould you like to remove it?")
            choice = input("y/n: ")            
            if str(choice).upper() == "Y":                
                cmd = 'rm -r '+folder
                subprocess.call(cmd, shell=True)
                cmd = 'mkdir '+folder
                subprocess.call(cmd, shell=True)                
                flag = False
                return True                
            elif str(choice).upper() == "N":                
                flag = False                
                return False                                  
            else:
                print("Please type y/Y or n/N")                               
    else:
        cmd = 'mkdir '+folder
        subprocess.call(cmd, shell=True)        
        return True
    

def tada_augmentation(X_train, y_train, obs_ids, i, tree_newick):
    """
    input: X_train, y_train, obs_ids, i, tree_newick
    """    
    X_train = np.transpose(X_train)    
    samp_ids = ['%d' % j for j in range(X_train.shape[1])]            
    otu_table = Table(X_train, obs_ids, samp_ids)
    ## Creates temp forlder to store each fold result
    folder = "tmp/"
    output = folder + "out"    
    if not os.path.exists(folder):        
        cmd = "mkdir -p {}".format(output)
        subprocess.check_output(cmd, shell = True)    
    
    fold_biom = folder + str(i) + ".biom"
    with biom.util.biom_open(fold_biom, "w") as f:
        otu_table.to_hdf5(f,"example")    
    labels_train = folder + "labels.tsv"
    pd.DataFrame(y_train).to_csv(labels_train, header = False, sep="\t")    
    
    # Execute TADA    
    cmd = "python TADA/src/utils/python/TADA_microbiom.py -t {} \
                        -b {} -o {} -g {}".format(tree_newick, fold_biom, output, labels_train)            
    subprocess.check_output(cmd, shell = True)    
    
    input_biom = output + "/" + str(i) + ".biom" 
    biom_tsv = output + "/" + str(i) + ".tsv"       
    cmd = "biom convert -i {} -o {} --to-tsv".format(input_biom, biom_tsv)
    subprocess.check_output(cmd, shell = True)    
    df_train_biom = pd.read_csv(biom_tsv, skiprows=[0], sep= "\t", index_col=0)
    
    input_aug_biom = output + "/augmented_data.biom"
    biom_aug_tsv = output + "/augmented_data.tsv"    
    cmd = "biom convert -i {} -o {} --to-tsv".format(input_aug_biom, biom_aug_tsv)
    subprocess.check_output(cmd, shell = True)    
    df_aug_biom = pd.read_csv(biom_aug_tsv, skiprows=[0], sep= "\t", index_col=0)
              
    labels = []
    tmp = [line.rstrip('\n') for line in open(output + "/labels.tsv")]
    for line in tmp:
        labels.append(int(line.split("\t")[-1]))
    tmp = [line.rstrip('\n') for line in open(output + "/augmented_meta_data.csv")]
    for line in tmp:
        labels.append(int(line.split("\t")[-1]))    
    # concatenate            
    X_train_tada = pd.concat([df_train_biom, df_aug_biom], ignore_index = True, axis=1).values
    y_train_tada = labels
    return np.transpose(X_train_tada), y_train_tada


def process_metrics_df(metrics_models, type_str):
    
    pd_type = pd.DataFrame.from_dict(metrics_models)
    pd_type = pd_type.stack().reset_index()
    tmp_mcc = pd_type[pd_type["level_1"].str.contains("mcc")]
    tmp_auc = pd_type[pd_type["level_1"].str.contains("auc")]
    tmp_auc["level_1"] = [i.split("_")[0] for i in tmp_auc.level_1]
    tmp_mcc["level_1"] = [i.split("_")[0] for i in tmp_mcc.level_1]
    pd_type = tmp_mcc
    pd_type.columns = ["fold", "model", "mcc"]
    pd_type["auc"] = tmp_auc[0].values
    pd_type["type"] = type_str    
    return  pd_type    



def smote_both(X_train, y_train):
                
    t = round(abs(min(dict(Counter(y_train)).values()) - max(dict(Counter(y_train)).values())))
    t_min = min(dict(Counter(y_train)).values())
    t_max = max(dict(Counter(y_train)).values())        
    over = SMOTE(sampling_strategy=(abs(t-t_min)/t_max))
    under = RandomUnderSampler(sampling_strategy=((t_max-t)/t_min))        
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)        
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote



def adasyn_both(X_train, y_train):
                
    t = round(abs(min(dict(Counter(y_train)).values()) - max(dict(Counter(y_train)).values())))    
    t_min = min(dict(Counter(y_train)).values())
    t_max = max(dict(Counter(y_train)).values())    
    over = ADASYN(sampling_strategy=(abs(t-t_min)/t_max))
    under = RandomUnderSampler(sampling_strategy=((t_max-t)/t_min))        
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)        
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote

