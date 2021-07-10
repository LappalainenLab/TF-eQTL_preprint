#!/bin/bash
##
##	ADASTRA_perms_wrapper.sh
##
## sbatch command:
##    num=$(wc -l TFs.info.curated.txt | awk '{print $1-1}')
##    job_name=runADASTRA_perms.R
##    sbatch --job-name=$job_name --array=1-$num --output=log/$job_name.%A_%a.out --open-mode=append --partition=pe2 --time=36:00:00 --mem=12G --mail-type=BEGIN,END,FAIL --mail-user=at3454@columbia.edu scripts/ADASTRA/ADASTRA_perms_wrapper.sh
##
##	ALT 01/19/2021
##

cd /gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB

# tf_file=TFs.info.curated.txt

# tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
#        tail -n 1 | awk '{print $1}'`

module load R/3.6.0
# Rscript scripts/ADASTRA/ADASTRA_perms.R $tf
Rscript scripts/ADASTRA/ADASTRA_perms.R

echo "Done."
date