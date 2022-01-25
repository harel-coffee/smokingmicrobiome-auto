
# Addressing class imbalance in microbiome-based classification using data augmentation and machine learning with smoking habit prediction from saliva as example


#### Celia Díez López, Diego Montiel González, Athina Vidaki, Manfred Kayser
Erasmus MC University Medical Center Rotterdam, Department of Genetic Identification, Rotterdam, The Netherlands.


## Requirements

#### Operating system: tested on Ubuntu 18.04LTS


### R: tested on R version 3.6.1 (2019-07-05) -- "Action of the toes"    
    Install dependencies, you can skip if already installed.
    
    gPCA v1.0
    dada2 v1.12.1
    ShortRead v1.42.0
    Biostrings v.2.52.0
    optparse v.16.4
    data.table v.1.12.6 
    tidyverse v.1.3.0
    rstatix v0.5.0
    ggpubr v0.3.0    
    DataCombine v0.2.21

### Python: Python >=3.6 with conda environment (recommended)
    
    $ pip install -r requirements.txt
    
    Installation using a conda environment (OPTIONAL)
    
    $ conda create --name <env> --file requirements.txt
    
    Activates environment
    $ conda activate <env>
    

### Tree-based Associative Data Augmentation (TADA) augmentation technique for classifying phenotypes based on the microbiome.
    
    $ git clone https://github.com/tada-alg/TADA

### External softwares (TADA preprocessing)
    MAFFT v7.310 
    FastTree v2.1.11


## Datasets
    
Datasets Study Accession numbers: PRJNA434300, PRJNA434312, PRJNA484874


* *data/eHOMD_v15.2_assignTaxonomy.fasta*: expanded Human Oral Microbiome Database (http://www.homd.org/)

* *data/metadata.csv*: Metadata table including phenotypes such as smoking status, sex and age.

* *cutadapt.R*: Remove adapters and primers for each of the two studies.

* *DADA2_pipeline*: Pipeline for taxaonomy assignation at species level for each study separately (https://benjjneb.github.io/dada2/tutorial.html).

* *batchEffect.R*: Assess for batch effects after dereplication and taxonomic assignation in Amplicon Sequence Variant (ASV) using guided Principal Component Analysis (gPCA).

### TADA preprocessing

* *ASV_consensus.R*: Generate consensus sequencing by multiple sequence alignemnt (DECIPHER) for each read taxonomic assignation (DADA2 pipeline).

* *TADA_preprocess.py*: Produce a rooted guided tree in newick format (FastTree).

* *NestedCV_pipeline.py*: Pipeline for Nested Cross-Validation using several Machine Learning models including: 
Logistic Regression, K-nearest neighbors, Decision Trees, Support Vector Machines Radial, Support Vector Machine Linear, Random Forest, XGBoosting
with additional data augmentation methods such as ADASYN, SMOTE, and TADA

* *stats.R*: Pairwise Wilcoxon test to compared data augmentation types.


contact at genid@erasmusmc.nl for any questions or issues concerning the scripts.




