
# Addressing class imbalance in microbiome-based classification using data augmentation and machine learning with smoking habit prediction from saliva as example


#### Celia Díez López, Diego Montiel González, Athina Vidaki, Manfred Kayser
Erasmus MC University Medical Center Rotterdam, Department of Genetic Identification, Rotterdam, The Netherlands.


## Requirements

#### Operating systems: Rested on Ubuntu 18.04LTS but can work on MacOS and Windows
#### Software: Make sure to install miniconda3 before to proceed. for additional information see the miniconda compatible to your system. https://docs.conda.io/en/latest/miniconda.html

    # Please run the following three steps (compulsory)
    # 1) Creates environment with external software, R and python packages
    conda create -n $PROJECT --file environment.txt -y
    # 2) Activate environment
    conda activate $PROJECT
    # 3) Install last R library
    Rscript -e 'devtools::install_github("christophergandrud/DataCombine")'

    # Deactivates environment (optional)
    conda deactivate
    

The study samples used can be found as the following Accession numbers: PRJNA434300, PRJNA434312, PRJNA484874

#### Additional data:

* *data/eHOMD_v15.2_assignTaxonomy.fasta*: expanded Human Oral Microbiome Database (http://www.homd.org/)

* *data/metadata.tsv*: Metadata table in TSV format including ID samples, and phenotypes such as smoking status, sex and age.

### Scripts and pipelines

* *cutadapt.R*: Remove adapters and primers for each of the two studies.

* *DADA2_pipeline*: Pipeline for taxaonomy assignation at species level for each study separately (https://benjjneb.github.io/dada2/tutorial.html).

* *batchEffect.R*: Assess for batch effects after dereplication and taxonomic assignation in Amplicon Sequence Variant (ASV) using guided Principal Component Analysis (gPCA).

#### TADA preprocessing

* *ASV_consensus.R*: Generate consensus sequencing by multiple sequence alignemnt (DECIPHER) for each read taxonomic assignation (DADA2 pipeline).

* *TADA_preprocess.py*: Produce a rooted guided tree in newick format (FastTree).


### Machine Learning Classifier workflow

* *NestedCV_pipeline.py*: Pipeline for Nested Cross-Validation using several Machine Learning models including: 
Logistic Regression, K-nearest neighbors, Decision Trees, Support Vector Machines Radial, Support Vector Machine Linear, Random Forest, XGBoosting
with additional data augmentation methods such as ADASYN, SMOTE, and TADA

* *stats.R*: Pairwise Wilcoxon test to compared data augmentation types.


contact at genid@erasmusmc.nl for any questions or issues concerning the scripts.




