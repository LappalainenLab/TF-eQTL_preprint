#!/bin/bash
##
##	calcaFCs.sh
##
##	EDF 4/15/2020
##

maf=$1 
tiss=$2
cov=$3
#tiss=`tail -n $SLURM_ARRAY_TASK_ID ../input_files/tissue_translation_colors_v8.txt | \
#	head -n 1 | awk -F'\t' '{print $2}'`	## note: tiss will be backwards
chr=$SLURM_ARRAY_TASK_ID

echo maf is $maf
echo tiss is $tiss
echo chr is chr$chr
if [[ "$cov" ]]; then
	echo aFC will be calculated with covariates
fi


phenotypes=phenotypes/filtered/by_chr/$tiss.v8.deseq_log2_expression.caviar_var.95set.eqtls.MAF$maf.overlap.chr$chr.bed.gz
genotypes=genotypes/by_chr/overlap_vars.MAF$maf.chr$chr.vcf.gz
if [[ "$cov" ]]; then
	covariates=covariates/GTEx_Analysis_v8_eQTL_covariates/$tiss.v8.covariates.txt
fi
qtls=../variant_sets/by_chr/caviar_var.95set.eqtls.MAF$maf.overlap.simple.chr$chr.txt

if [[ "$cov" ]]; then
	output_file=by_tiss/$tiss/$tiss.afcs.cov.overlap_vars.MAF$maf.chr$chr.txt
else
	output_file=by_tiss/$tiss/$tiss.afcs.nocov.overlap_vars.MAF$maf.chr$chr.txt
fi

echo output file will be $output_file

source ~/.bashrc

aFC.py --vcf $genotypes \
	--pheno $phenotypes \
	--qtl $qtls \
	`if [[ "$cov" ]]; then echo --cov $covariates; fi` \
	--log_xform 1 --log_base 2 \
	--o $output_file

wc -l $qtls
wc -l $output_file

echo "Done."
date




