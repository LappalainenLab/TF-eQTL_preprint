#!/bin/bash
##
##	formatTFCorrsMed0_byTF_wrapper.sh
##
##	Script to run formatTFCorrsMed0_byTF.R
##
##	EDF 6/17/20
##

cd ~/projects/crosstiss_tf_corrs/

tf_file=input_files/TFs.info.curated.txt

tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
        tail -n 1 | awk '{print $1}'`

module load R/3.6.0
Rscript scripts/formatTFCorrsMed0_byTF.R $tf

echo "Done."
date



