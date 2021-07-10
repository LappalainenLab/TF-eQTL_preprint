#!/bin/bash
##
##	Submits getExpr script 
##		for each tissue
##
##	EDF 4/13/2020
##


maf=$1
tiss=`tail -n $SLURM_ARRAY_TASK_ID ../input_files/tissue_translation_colors_v8.txt | \
	head -n 1 | awk -F'\t' '{print $2}'`    ## note: tiss will be backwards
echo $maf
echo $tiss

mkdir by_tiss/$tiss

job_name=calcaFC_cov_$tiss
sbatch --job-name=$job_name --array 1-22 \
        --output=log/$job_name.%A.%a.out \
        --partition=pe2 --time=100:00:00 \
        --mail-type=BEGIN,END,FAIL --mail-user=eflynn@nygenome.org \
        scripts/calcaFC.sh $maf $tiss cov

job_name=calcaFC_nocov_$tiss
sbatch --job-name=$job_name --array 1-22 \
        --output=log/$job_name.%A.%a.out \
        --partition=pe2 --time=100:00:00 \
        --mail-type=BEGIN,END,FAIL --mail-user=eflynn@nygenome.org \
        scripts/calcaFC.sh $maf $tiss




