#!/bin/bash
##
##	permTFEnrich_bytf_wrapper.sh
##
##
##	EDF 6/22/2020
##


cd ~/projects/overlap/

gene_list=$1
tf_file=input_files/TFs.info.curated.txt

tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
        tail -n 1 | awk '{print $1}'`

echo $gene_list $tf

module load R/3.6.0
Rscript scripts/permGenesTFs_bytf3.R $gene_list $tf

echo "Done."
date







