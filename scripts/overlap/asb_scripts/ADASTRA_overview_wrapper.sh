#!/bin/bash
##
##	ADASTRA_overview_wrapper.sh
##
##	Script to run ADASTRA_overview.R
##
## sbatch command:
##    num=$(wc -l TFs.info.curated.txt | awk '{print $1-1}')
##    job_name=runADASTRA_overview
##    sbatch --job-name=$job_name --array=1-$num --output=log/$job_name.out --open-mode=append --partition=pe2 --time=100:00:00 --mem=12G --mail-type=BEGIN,END,FAIL --mail-user=at3454@columbia.edu scripts/ADASTRA/ADASTRA_overview_wrapper.sh

cd /gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB

tf_file=TFs.info.curated.txt

tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
        tail -n 1 | awk '{print $1}'`

module load R/3.6.0
Rscript scripts/ADASTRA/ADASTRA_overview.R $tf $fdr

echo "Done."
date