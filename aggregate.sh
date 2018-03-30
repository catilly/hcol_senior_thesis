#!/bin/bash

######

METHOD=$1 # One of marginal, PRS, BLUP, lasso, EN, ENcombined, XPBLUP, or XPEN; it should match the run.METHOD_... filenames
DIR=$2 # Full path to directory containing the output files of run_all.sh; these files will be combined

######
# Example: sh aggregate.sh XPBLUP /bcb/agusevlab/cali/GEUVADIS/Exp49_20180216/GeneQuant
######

# Get rid of any empty files so that they don't mess with the headers
find ${DIR}/run.METHOD_${METHOD}.*.out -size 0 -delete

# Concatenate files
# https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
if [ -f ${DIR}/concat.METHOD_${METHOD}.out ]; then
    echo "File already exists."
else
    array=( ${DIR}/run.METHOD_${METHOD}.*.out )
    head -1 ${array[0]} > ${DIR}/concat.METHOD_${METHOD}.out 
    tail -n +2 -q ${array[@]:0} >> ${DIR}/concat.METHOD_${METHOD}.out
    rm ${DIR}/run.METHOD_${METHOD}.*.out
fi

# Process for plotting
/apps/R-3.4.1share/bin/R --slave --args "${DIR}/concat.METHOD_${METHOD}.out" $METHOD < /bcb/agusevlab/cali/GEUVADIS/process_files_for_plotting_GEUVADIS.R > ${DIR}/long.METHOD_${METHOD}.out