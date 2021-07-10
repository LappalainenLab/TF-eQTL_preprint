#!/bin/bash
##
##	getCaviarMAF.sh
##
##	Get variants of a certain MAF cutoff
##		that were in a Caviar 95% conf set.
##
##	EDF 3/30/20
##

num=00`echo $SLURM_ARRAY_TASK_ID-1 | bc`
num_str=${num: -3}
cav_bed=caviar_output_GTEx_LD_combined/all_vars_split/split.$num_str
#cav_bed=caviar_output_GTEx_LD_combined/all_tiss.all_snps.95set.bed
maf=$1
maf_vcf=caviar_output_maf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF$maf.noX.vcf.gz
out_vcf=caviar_output_maf/split/caviar_var.95set.all_tiss.MAF$maf.$num_str.vcf

#tabix -H $maf_vcf > $out_vcf
while read line
do
	`echo $line | awk -v file=$maf_vcf \
		'{print "tabix",file,$1":"$3"-"$3}'`
done < $cav_bed | tee $out_vcf.temp | \
	sort | uniq > $out_vcf
wc -l $out_vcf
srm $out_vcf.temp
#bgzip -f $out_vcf
#tabix -f $out_vcf.gz



