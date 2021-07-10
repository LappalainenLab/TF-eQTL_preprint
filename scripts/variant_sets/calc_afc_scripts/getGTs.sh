#!/bin/bash
##
##	getGTs.sh
##
##	EDF 4/15/20
##


chri=$SLURM_ARRAY_TASK_ID
variant_set=$1

echo variant set is $variant_set
echo chr is $chri

input_vcf=genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz
output_vcf=genotypes/by_chr/$variant_set.chr$chri.vcf

echo output vcf will be $output_vcf

gunzip -c $input_vcf | head -n 5000 | grep "CHROM" | head -n 1 > $output_vcf
while read line
do
	snp=`echo $line | awk '{print $3}'`
	#echo $snp
        chr=`echo $line | awk '{print $1}'`
        pos=`echo $line | awk '{print $2}'`
        tabix $input_vcf $chr\:$pos\-$pos | grep $snp
done < ../variant_sets/by_chr/$variant_set.chr$chri.txt >> $output_vcf

wc -l ../variant_sets/by_chr/$variant_set.chr$chri.txt
wc -l $output_vcf

bgzip -f $output_vcf
tabix -f $output_vcf.gz

echo "Done."
date



