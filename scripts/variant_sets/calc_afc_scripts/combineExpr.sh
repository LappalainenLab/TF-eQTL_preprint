#!/bin/bash
##
##	Pulls genes of interest from 
##		expression file for
##		the designated tissue
##
##	EDF 4/13/2020
##


variant_set=$1 ## base stem without .txt or ending
tiss=`tail -n $SLURM_ARRAY_TASK_ID ../input_files/tissue_translation_colors_v8.txt | \
        head -n 1 | awk -F'\t' '{print $2}'`    ## note: tiss will be backwards
echo $variant_set
echo $tiss

orig_file=phenotypes/afc_deseq_log2_expression_allgenes/$tiss.v8.deseq_log2_expression.bed.gz
out_file=phenotypes/filtered/$tiss.v8.deseq_log2_expression.$variant_set.bed
zcat $orig_file | head -n 1 > $out_file
for chr in `seq 1 22`
do
	in_file_chr=phenotypes/filtered/by_chr/$tiss.v8.deseq_log2_expression.$variant_set.chr$chr.bed 
	cat $in_file_chr >> $out_file
done

bgzip -f $out_file
tabix -f -p bed $out_file.gz

echo "Done."
date




