# Supplementary code for "Novel methods for cross-population prediction of gene expression"
Last updated: April 1, 2018, by Catherine Li
This repository contains supplementary code for my honors undergraduate thesis, which was advised by Alexander Gusev and submitted to the Harvard Statistics Department. Here, we set up a walkthrough for evaluating eight methods -- top SNP, polygenic risk score, BLUP, lasso, EN, combined EN, XP-BLUP, and XP-EN -- on 1000 Genomes/Geuvadis data from chromosome 21. Beyond illustrating our workflow, we document a suite of custom R functions and scripts here that enable data processing, prediction, and simulation.

## Description of Contents

This directory contains the following sub-directories and files:

- **Data**: Contains preprocessed data for chromosome 21.
    - **Data/C1_lists**: Contains SNP lists for chr21 genes for use in XP-BLUP and XP-EN. These were obtained by running association tests in the EUR data using the plink --assoc command and a p-value cutoff of 1.5e-4.
    - **Data/Expression**: Contains Geuvadis expression data from Lappalainen et al. (2013). The data has been filtered to include only chr21 genes for which Lappalainen et al. discovered a significant association in either EUR or YRI.
    - **Data/Genotype**: Contains 1000 Genomes ch21 genotype data. Duplicate SNPs have been removed. 1000G.POP records the population of each individual in the dataset.

- **out**: An empty directory where output files will be stored.

- **XP-BLUP_v0.3**: Contains files and scripts released by Coram et al. (2017) that implement XP-BLUP. Aside from **XP-BLUP_v0.3/xpblup-CL.sh**, which is a slightly modified version of the original release that makes parallelization easier, all of these files can be obtained from http://med.stanford.edu/tanglab/software/XPBLUP/XPBLUP.html.

- **aggregate.sh**: A script for aggregating prediction scores across genes.

- **predict_combined.R**: A script for fitting models on a combined EUR and YRI training set and then predicting in a subset of the YRI population. This is used to implement combined EN.

- **predict_EUR_only_methods.R**: A script for fitting models in EUR and predicting in YRI. This is used to implement top SNP, PRS, BLUP, lasso, and EN.

- **predict_XP-EN.R**: A script for fitting and evaluating XP-EN in YRI.

- **process_files_for_plotting.R**: Used in aggregate.sh to process files into a format that can easily be visualized in ggplot2.

- **run_all.sh**: A script that runs all eight methods on the chr21 data.

- **run_XP-BLUP.sh**: A script for fitting and evaluating XP-BLUP in YRI.

- **simulate_expression_prediction_v2.R**: A set of custom R functions that are used throughout our study, documented in detail. This contains processing and model-fitting functions (including our novel method XP-EN) employed in the above scripts, but it also contains functions for simulating gene expression.

- **XP-BLUP_helper.R**: A script that generates appropriately formatted input files for XP-BLUP.

## Demo: Evaluating Expression Prediction Methods in Chromosome 21

In this demo, we will run and evaluate eight methods -- top SNP (also called marginal), polygenic risk score (PRS), BLUP, lasso, EN, combined EN, XP-BLUP, and XP-EN -- on genes in chromosome 21.

### Prerequisites

TODO

Before we begin, note that this demo requires

- **R** v.3.4.1 or above (https://www.r-project.org/), with the following packages installed:
    - snpStats
    - glmnet
    - plyr
    - dplyr
    - reshape2
    - MASS

- **plink**

- **gcta**

### Running the demo

After downloading all of the requisite software, move into this directory and run
```
sh run_all.sh [R_LOC] [PLINK_LOC] [GCTA_LOC]
```
from the Unix command line. R_LOC, PLINK_LOC, and GCTA_LOC should be the paths to R, plink, and gcta respectively in your system. This script will run all eight methods and output files into the out/ directory. It can take a while, so you may want to split up the script and parallelize it on a computing cluster, or simply comment out the methods that are not of interest.

Once this script is finished running, the output files can be concatenated and processed by typing
```
sh aggregate.sh [METHOD] [DIR]
```
where DIR is the location of the files you plan to aggregate and METHOD is method for which you are aggregating prediction scores (one of marginal, PRS, BLUP, lasso, EN, XPBLUP, or XPEN -- this argument should match the naming scheme of the files in DIR).

You resulting files, e.g. "long.METHOD_XPBLUP.out," contain the predictive correlation achieved by a particular method on each gene in chromosome 21. They can be loaded into the data analysis software of your choice and further explored.
