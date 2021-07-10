#!/bin/bash
##
##	submitCorrelateTFs_med0_wrapper.sh
##
##	Script to run correlateTFs.R
##
##	EDF 4/29/20
##

cd ~/projects/crosstiss_tf_corrs/

tf_file=input_files/TFs.info.curated.txt

maf=$1
tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
	tail -n 1 | awk '{print $1}'`

#module load R/3.6.0
#Rscript scripts/correlateTFs.R $maf $tf
job_name=crossTissCorrs_med0_$tf\_05
sbatch --job-name=$job_name --array 1-22 \
	--output=log/$job_name.%A.%a.out \
	--partition=pe2 --time=100:00:00 \
	--mem=8G \
	--mail-type=BEGIN,END,FAIL --mail-user=eflynn@nygenome.org \
	scripts/correlateTFs_med0_bychr_wrapper.sh $maf $tf

echo "Done."
date


