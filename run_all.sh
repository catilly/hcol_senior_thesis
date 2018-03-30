#!/bin/bash

#####

R_LOC=$1 # Location of R on your system (https://www.r-project.org, v3.4.1 or later)
PLINK_LOC=$2 # Location of plink (https://www.cog-genomics.org/plink2, v1.90 or later)
GCTA_LOC=$3 # Location of gcta (http://cnsgenomics.com/software/gcta, v1.25.3 or later)

QUANT_DATA=Data/Expression/GeneQuantRPKM.filt.chr21.txt # Full path to quantification data table, assumed to be from a particular chromosome
ORIG_PLINK=Data/Genotype/1000G.all.deduped.21 # Full path to plink genotype data; should be from the same chromosome as QUANT_DATA
C1_LOC=Data/C1_lists # Directory where C1 lists are stored
OUT_PATH=out # Full path to directory where output will be stored

######
# Example: sh run_all.sh /apps/R-3.4.1share/bin/R /bcb/agusevlab/cali/plink_linux_x86_64/plink /bcb/agusevlab/cali/gcta_1.26.0/gcta64
######

# Get chromosome number from input file
CHR=$(echo $QUANT_DATA| cut -d'.' -f 3)

# Run top SNP, PRS, BLUP, lasso, and EN
for m in marginal PRS BLUP lasso EN; do
  $R_LOC --slave --args $QUANT_DATA $m < predict_EUR_only_methods.R > $OUT_PATH/run.METHOD_$m.${CHR}.out;
done

# Run combined EN
$R_LOC --slave --args $QUANT_DATA EN < predict_combined.R > $OUT_PATH/run.METHOD_ENcombined.${CHR}.out

# Run XP-EN
$R_LOC --slave --args $QUANT_DATA $C1_LOC < predict_XP-EN.R > $OUT_PATH/run.METHOD_XPEN.${CHR}.out

# Run XP-BLUP
sh run_XP-BLUP.sh $R_LOC $PLINK_LOC $GCTA_LOC $QUANT_DATA $ORIG_PLINK $C1_LOC XP-BLUP_helper.R $OUT_PATH