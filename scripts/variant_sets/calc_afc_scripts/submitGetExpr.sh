#!/bin/bash
##
##	Submits getExpr script 
##		for each tissue
##
##	EDF 4/13/2020
##


variant_set=$1 ## base stem without .txt or ending
tiss=`tail -n $SLURM_ARRAY_TASK_ID ../input_files/tissue_translation_colors_v8.txt | \
	head -n 1 | awk -F'\t' '{print $2}'`    ## note: tiss will be backwards
echo $variant_set
echo $tiss

job_name=getExpr_$tiss
sbatch --job-name=$job_name --array 1-22 \
        --output=log/$job_name.%A.%a.out \
        --partition=pe2 --time=100:00:00 \
        --mail-type=BEGIN,END,FAIL --mail-user=eflynn@nygenome.org \
        scripts/getExpr.sh $variant_set $tiss


