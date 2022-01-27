
# Addressing class imbalance in microbiome-based classification using data augmentation and machine learning with smoking habit prediction from saliva as example


#### Celia Díez López, Diego Montiel González, Athina Vidaki, Manfred Kayser
Erasmus MC University Medical Center Rotterdam, Department of Genetic Identification, Rotterdam, The Netherlands.


    This repository contains a general Machine Learning Classifier workflow meant to be used for microbiome data, where the sample size (n) of the phenotypes/labels/classes of interest is unbalanaced. 
    You can select among several Machine Learning Classifiers to be used in your data. 
    We implemented a Nested Cross-Validation method together with augmentation techniques such as SMOTE, ADASYN, and [TADA](https://github.com/tada-alg/TADA/) [https://doi.org/10.1093/bioinformatics/btz394]


## Requirements

#### Operating systems: Tested on Ubuntu 18.04LTS but can work on MacOS and Windows
#### Software: Make sure to install miniconda3 before to proceed. for additional information see the miniconda compatible to your system. https://docs.conda.io/en/latest/miniconda.html

    # Please run the following three steps (compulsory)
    # 1) Creates environment with external software, R and python packages
    conda create -n $PROJECT --file environment.txt -y
    # 2) Activate environment
    conda activate $PROJECT
    
    # Deactivates environment (optional)
    conda deactivate
    

The study samples used can be found as the following Accession numbers: PRJNA434300, PRJNA434312, PRJNA484874

#### Additional data:

* *data/eHOMD_v15.2_assignTaxonomy.fasta*: expanded Human Oral Microbiome Database (http://www.homd.org/)

* *data/metadata.tsv*: Metadata table in TSV format including ID samples, and phenotypes such as smoking status, sex and age.

### Scripts and pipelines

#### TADA preprocessing (optional)

    TADA is a 

* *ASV_consensus.R*: Generate consensus sequencing by multiple sequence alignemnt (DECIPHER) for each read taxonomic assignation (DADA2 pipeline).

* *TADA_preprocess.py*: Produce a rooted guided tree in newick format (FastTree).


#### Machine Learning Classifier workflow

* *NestedCV_pipeline.py*: Pipeline for Nested Cross-Validation using several Machine Learning models including: 
Logistic Regression, K-nearest neighbors, Decision Trees, Support Vector Machines Radial, Support Vector Machine Linear, Random Forest, XGBoosting
with additional data augmentation methods such as ADASYN, SMOTE, and TADA

* *stats.R*: Pairwise Wilcoxon test to compared data augmentation types.


contact at genid@erasmusmc.nl for any questions or issues concerning the scripts.




