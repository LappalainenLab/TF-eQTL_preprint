#!/bin/bash
##
##
##	EDF 2/25/2021
##


cd ~/projects/IRF1/

module load tabix


out_vars=coding_snps/all_gene_exon_var_gt.txt
out_vars_het=coding_snps/all_gene_exon_var_gt.hets.txt
echo gene exon variant chr pos ref alt gt | tr ' ' '\t' | \
        tee $out_vars > $out_vars_het


for file in coding_snps/exons/*exons.txt
do
        gene=`basename $file | awk -F. '{print $1"."$2}'`
	echo $gene 1>&2

        while read line
        do
                chr=`echo $line | awk '{print $1}'`
                pos1=`echo $line | awk '{print $4}'`
                pos2=`echo $line | awk '{print $5}'`
                exon=`echo $line | tr ';' '\n' | grep exon_number | awk '{print $2}' | sed 's/"//'`
                gt_line=`tabix input_files/293_CG.liftoverhg38.sorted.cleaned.vcf.gz $chr":"$pos1"-"$pos2`

        if [[ "`echo $gt_line`" ]]
        then
                var_id=`echo $gt_line | awk '{print $1"_"$2"_"$4"_"$5"_b38"}'`
                gt=`echo $gt_line | awk '{print $10}'`
                echo $gene $exon $var_id `echo $var_id | tr '_' ' ' | sed 's/b38//'`$gt | \
                         tr ' ' '\t'
        fi

        done < $file
done | tee -a $out_vars | grep "0/1" >> $out_vars_het

cat $out_vars_het | awk '{print $4,$5}' > coding_snps/all_coding_het_pos.txt




