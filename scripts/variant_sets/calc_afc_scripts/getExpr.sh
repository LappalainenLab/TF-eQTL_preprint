#!/bin/bash
##
##	Pulls genes of interest from 
##		expression file for
##		the designated tissue
##
##	EDF 4/13/2020
##

variant_set=$1 ## base stem without .txt or ending
tiss=$2
#tiss=`tail -n $SLURM_ARRAY_TASK_ID ../input_files/tissue_translation_colors_v8.txt | \
#	head -n 1 | awk -F'\t' '{print $2}'`	## note: tiss will be backwards
chr=$SLURM_ARRAY_TASK_ID

echo variant set is $variant_set
echo tiss is $tiss
echo chr is chr$chr
in_expr_file=phenotypes/afc_deseq_log2_expression_allgenes/$tiss.v8.deseq_log2_expression.bed.gz
out_expr_file=phenotypes/filtered/by_chr/$tiss.v8.deseq_log2_expression.$variant_set.chr$chr.bed

gene_grep=`awk '{print $1}' ../variant_sets/$variant_set.genesort.chr$chr.txt | \
	uniq | tr '\n' '|' | sed 's/[|]$/|GTEX\n/'`

zgrep -E $gene_grep $in_expr_file | sort -k 1,1 -k 2n,2n -k 3n,3n  > $out_expr_file
#bgzip -f $out_expr_file
#tabix -f -p bed $out_expr_file.gz

echo "Done."
date




