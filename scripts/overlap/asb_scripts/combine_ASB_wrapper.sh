#!/bin/bash
##
##	ADASTRA_perms_wrapper.sh
##
## sbatch command:
##    num=60
##    job_name=runCombine_ASB.R
##    sbatch --job-name=$job_name --array=0-$num --output=log/$job_name.out --open-mode=append --partition=pe2 --time=10:00:00 --mem=12G --mail-type=BEGIN,END,FAIL --mail-user=at3454@columbia.edu scripts/ADASTRA/combine_ASB_wrapper.sh
##
##	ALT 01/19/2021
##

cd /gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB

module load R/3.6.0
Rscript scripts/ADASTRA/combine_ASB.R $SLURM_ARRAY_TASK_ID

echo "Done."
date