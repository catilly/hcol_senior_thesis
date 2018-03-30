#!/bin/bash

#####

R_LOC=$1 # Location of R on your system (https://www.r-project.org, v3.4.1 or later)
PLINK_LOC=$2 # Location of plink (https://www.cog-genomics.org/plink2, v1.90 or later)
GCTA_LOC=$3 # Location of gcta (http://cnsgenomics.com/software/gcta, v1.25.3 or later)

QUANT_DATA=$4 # Full path to quantification data table, assumed to be from a particular chromosome
ORIG_PLINK=$5 # Full path to plink genotype data; should be from the same chromosome as QUANT_DATA
C1_LOC=$6 # Directory where C1 lists are stored
IN_PATH=$7 # # Full path to R script for preprocessing
OUT_PATH=$8 # Full path to directory where output will be stored

#####

# Create SNP lists and phenotype files for pop 1, pop 2 train, and pop 2 test
$R_LOC --slave --args $QUANT_DATA $OUT_PATH < $IN_PATH

# Get chromosome number from input file
CHR=$(echo $QUANT_DATA| cut -d'.' -f 3)

# Get list of TargetIDs from input file
ID_LIST=$(awk '{print $1}' ${QUANT_DATA})

# Loop over genes
for i in ${ID_LIST}; do

  # Create relevant plinks
  for PHENO in ${OUT_PATH}/pop1.pheno.${CHR}.TargetID_${i} ${OUT_PATH}/pop2.train.pheno.${CHR}.TargetID_${i} ${OUT_PATH}/pop2.test.pheno.${CHR}.TargetID_${i}; do
    $PLINK_LOC --bfile $ORIG_PLINK --extract ${OUT_PATH}/locus.SNP.list.${CHR}.TargetID_${i}.txt --keep ${PHENO}.txt --pheno ${PHENO}.txt --make-bed --allow-no-sex --out ${PHENO}
  done

#  # If not already done, create C1 lists by running an association test on the pop1 file
#  # Note that the p-value threshold is set to 1.5e-4, which is approximately the Bonferroni correction for the average number of SNPs in 1 MB of chr1
#  /bcb/agusevlab/cali/plink_linux_x86_64/plink --bfile ${OUT_PATH}/pop1.pheno.${CHR}.TargetID_${i} --assoc --allow-no-sex --pfilter 0.00015 --out ${OUT_PATH}/assoc.${CHR}.TargetID_${i}
#  tail -n +2 ${OUT_PATH}/assoc.${CHR}.TargetID_${i}.qassoc | sed s/^\ *//g | tr -s ' ' ' ' | cut -f 2 -d ' ' > ${OUT_PATH}/C1.${CHR}.TargetID_${i}.txt
#  rm ${OUT_PATH}/pop1.pheno.${CHR}.TargetID_${i}*
#  rm ${OUT_PATH}/assoc.${CHR}.TargetID_${i}*

  # XP-BLUP
  # All causal variants
  XP-BLUP_v0.3/xpblup-CL.sh --plink=$PLINK_LOC --gcta=$GCTA_LOC --train=${OUT_PATH}/pop2.train.pheno.${CHR}.TargetID_${i} --test=${OUT_PATH}/pop2.test.pheno.${CHR}.TargetID_${i} --snplist=${C1_LOC}/C1.${CHR}.TargetID_${i}.txt --outdir=${OUT_PATH} --outprefix=out.${CHR}.TargetID_${i}
  # Clean up
  awk -v d="${i}" 'BEGIN{FS=OFS="\t"} {print $0, (NR>1?d:"CORR")}' ${OUT_PATH}/out.${CHR}.TargetID_${i}.predict.profile > ${OUT_PATH}/final_preds.${CHR}.TargetID_${i}.txt

done

# More cleanup
rm ${OUT_PATH}/*pheno.${CHR}.TargetID_*
# rm ${OUT_PATH}/C1.${CHR}.TargetID_*
rm ${OUT_PATH}/out.${CHR}.TargetID_*
rm ${OUT_PATH}/locus.SNP.list.${CHR}.TargetID_*
find ${OUT_PATH}/final_preds.${CHR}.TargetID_* -size 0 -delete # Delete empty files; sometimes the REML fails

# Aggregate across correlations
# https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
if [ -f ${OUT_PATH}/run.METHOD_XPBLUP.${CHR}.out ]; then
  echo "File already exists."
else
  array=( ${OUT_PATH}/final_preds.${CHR}.TargetID_* )
  head -1 ${array[0]} > ${OUT_PATH}/run.METHOD_XPBLUP.${CHR}.out
  tail -n +2 -q ${array[@]:0} >> ${OUT_PATH}/run.METHOD_XPBLUP.${CHR}.out
  rm ${OUT_PATH}/final_preds.${CHR}.TargetID_*
fi