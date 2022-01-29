
# Addressing class imbalance in microbiome-based classification using data augmentation and machine learning with smoking habit prediction from saliva as example


#### Celia Díez López, Diego Montiel González, Athina Vidaki, Manfred Kayser
Erasmus MC University Medical Center Rotterdam, Department of Genetic Identification, Rotterdam, The Netherlands.

> This method helps you to evalute your microbiome data, and further investigate if your phenotype of interest has the prediction potential to be used in a Classifier.


This repository contains a general Machine Learning Classifier workflow meant to be used for microbiome data, where the sample size *(n)* of the phenotypes/labels/classes of interest is unbalanaced. 
You can select among several Machine Learning Classifiers to be used in your data. 
We implemented a Nested Cross-Validation method together with augmentation techniques such as SMOTE, ADASYN,and [TADA](https://github.com/tada-alg/TADA/) [https://doi.org/10.1093/bioinformatics/btz394]


## Requirements

#### Operating systems: Tested on Ubuntu 18.04LTS but can work on MacOS and Windows
#### Software: Make sure to install miniconda3 before to proceed. for additional information see the miniconda compatible to your system. https://docs.conda.io/en/latest/miniconda.html

    # Please run the following three steps (compulsory)
    # 1) Creates environment with external software, R and python packages
    $ conda create -n $PROJECT --file environment.txt -y
    # 2) Activate environment
    $ conda activate $PROJECT
    
    # Deactivates environment (optional)
    $ conda deactivate
    

The study samples used can be found as the following Accession numbers: PRJNA434300, PRJNA434312, PRJNA484874

#### Additional data:

* *data/eHOMD_v15.2_assignTaxonomy.fasta*: expanded Human Oral Microbiome Database (http://www.homd.org/)

* *data/metadata.tsv*: Metadata table in TSV format including ID samples, and phenotypes such as smoking status, sex and age.

## Scripts and pipelines

In order to proceed you can use the data provided in the repository or use your own. For your own data we recommend an input provided by [DADA2](https://benjjneb.github.io/dada2/tutorial.html). We recommend to follow their tutorial to obtain the Amplicon Sequence Variant (ASV) file (e.g. Genus or Species) and The Sequence from each of the taxas assignated to each ASV.

## TADA preprocessing (optional)

These following steps are optional, are needed to run augmentation with TADA in NestedCV_pipeline.py. 

### **ASV_consensus.R: Generate consensus sequencing by multiple sequence alignemnt (DECIPHER) for each read taxonomic assignation (DADA2 pipeline).**

Help function
```
$ Rscript --vanilla ASV_consensus.R -h

Usage: ASV_consensus.R [options]


Options:
	-s CHARACTER, --species=CHARACTER
		ASV at the Species-level, output from DADA2

	-t CHARACTER, --taxas=CHARACTER
		Input file of assigned taxas from DADA2 pipeline, at the
              Species-level

	-o CHARACTER, --output=CHARACTER
		output folder

	-h, --help
		Show this help message and exit

```
*Example*

```
$ Rscript --vanilla ASV_consensus.R -s data/species.csv -t data/taxas_eHMOD.csv -o output
```
The following file is generated:

* $output/Species_OTU_seqs_consensus.csv: Consensus of taxaonomic sequences using the [_ConsensusSequence_] (https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html) function from DECIPHER Bioconductor package.


### **TADA_preprocess.py: Produce a rooted guided tree in newick format (FastTree).**


Help function
```
$python TADA_preprocess.py -h

usage: TADA_preprocess.py [-h] -c FILE -s FILE -o DIR

optional arguments:
  -h, --help            show this help message and exit
  -c FILE, --consensus FILE
                        output/Species_OTU_seqs_consensus.csv
  -s FILE, --species FILE
                        data/species.csv
  -o DIR, --tada_out DIR
                        output/TADA/
```


*Example*

```
$ python TADA_preprocess.py -c output_3/Species_OTU_seqs_consensus.csv -s data/species.csv -o output_3/TADA
```

The following files are generated:

* $output/TADA/taxa_species.fasta: Taxonomic sequences in FASTA format.
* $output/TADA/taxa_species_mafft.fasta: Taxonomic sequences after Multiple Sequence Alignment with MAFFT v7.490.
* $output/TADA/taxa_species_tree.nwk: Approximately-Maximum-Likelihood constructed tree based on the align sequences using FastTree version 2.1.10.




### **Machine Learning Classifier workflow**

* *NestedCV_pipeline.py*: Pipeline for Nested Cross-Validation using several Machine Learning models including: 
Logistic Regression, K-nearest neighbors, Decision Trees, Support Vector Machines Radial, Support Vector Machine Linear, Random Forest, XGBoosting
with additional data augmentation methods such as ADASYN, SMOTE, and TADA


Here is the help function:
```
$ python NestedCV_pipeline.py -h

usage: NestedCV_pipeline.py [-h] -i FILE -o PATH -m FILE -target STRING -id
                            STRING [-t INT] [-r FILE] [-ml MODEL]
                            [-a AUGMENTED] [-iter INT]

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        table taxa species e.g. data/species_intersect.csv.
                        Taxa as columns, sample names as rows
  -o PATH, --output PATH
                        output folder name e.g. data/output/
  -m FILE, --metadata FILE
                        table with metadata data/metadata.csv
  -target STRING, --target STRING
                        Target phenotype name based on the column of metadata
                        e.g. smoking
  -id STRING, --id STRING
                        Column id sample names on the column of metadata e.g.
                        SRR or Samples
  -t INT, --threads INT
                        threads to be used during Nested Cross-Validation e.g.
                        4
  -r FILE, --root-tree FILE
                        Rooted tree in newick format from ASV consensus e.g.
                        data/TADA/taxa_species_tree.nwk
  -ml MODEL, --ml-model MODEL
                        Machine learning models, if more than one is selected
                        implements NestedCV and reports output and prediction
                        estimations using MCC and AUC, otherwise you can
                        select one to export as final model (also needs to be
                        selected an augmented data type) [LR, KNN, DT, SVMR,
                        SVML, RF, XG] e.g. LR,KNN,RF or e.g. SVML (will
                        generate one model and export)
  -a AUGMENTED, --augmented AUGMENTED
                        Data type augmented technique, if no type is selected
                        uses all data types [DEFAULT, ADASYN_over,
                        ADASYN_both, SMOTE_both, SMOTE_over, TADA]
  -iter INT, --iterations INT
                        number of iterations [DEFAULT:10]
```

*Examples*

```
# 1) To run with all models and augmentation techniques:
$ python NestedCV_pipeline.py -i data/species.csv -o output_3 -m data/metadata.tsv -target smoking_status -id SRR -t 4 -r output/TADA/taxa_species_tree.nwk -iter 10

# 2) To run with two models and two augmentation techniques:
$ python NestedCV_pipeline.py -i data/species.csv -o output_3 -m data/metadata.tsv -target smoking_status -id SRR -t 4 -r output/TADA/taxa_species_tree.nwk -iter 10
```
From step 1 and 2 the following files are generated:

* $output/mean_results_train.txt: Table of mean probabilities used on each combination of augmentation and machine learning model technique, AUC and MCC with _std_ (standard deviation)
* $output/mcc_boxplot.png: Boxplot containing same values of mean_results_train.txt but with Matthews Correlation Coefficient.
* $output/auc_boxplot.png: Boxplot containing same values of mean_results_train.txt but with The Aurea Under the Curve.

#### Note: Step 1 or 2 allows you to visualize prediction score using Matthews Correlation Coefficient (MCC) and The Area-Under the Curve (AUC). Once you have chosen the best option of model + augmentation technique you can use it in step 3 to train only with your best combination and see predicted values with a random subset (20%) of your data (unseen). For more information please see scheme from the paper about Nested Cross-validation approach.  

```
# 3) Select one model and augmentation technique.
$ python NestedCV_pipeline.py -i data/species.csv -o output -m data/metadata.tsv -target smoking_status -id SRR -t 7 -r output_3/TADA/taxa_species_tree.nwk -ml KNN -a TADA -iter 10
```
From step 3 the following files are generated:

* $output/MODEL/: This folder will contain the models use from step 3. Given the augmentation and machine learning model technique. 
* $output/MODEL/probs_test.txt: Predicted probabilities used in 20% of random samples from your data used as test against the rest of the models used in the Nest-Cross validation approach.



contact at genid@erasmusmc.nl for any questions or issues concerning the scripts.




