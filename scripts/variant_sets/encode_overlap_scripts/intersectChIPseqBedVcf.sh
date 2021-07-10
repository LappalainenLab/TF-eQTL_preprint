#!/bin/bash

source ~/.bashrc
module load bedtools
module load tabix
cd /commons/groups/lappalainen_lab/data/gtex/v8/encode_overlap/

in_vcf=$1
in_bed=$2
in_tf=$3
in_cell=$4
out_vcf=peak_intersect_vcfs/$in_tf.$in_cell.`basename $in_bed | sed 's/bed.gz//' | sed 's/sorted.merged.//'``basename $in_vcf | sed 's/.gz//'`


echo VCF is $in_vcf
echo bed file is $in_bed
echo tf is $in_tf
echo cell type is $in_cell
echo output file will be written to $out_vcf.gz

bedtools intersect -a $in_vcf -b $in_bed -wa | cut -f 1-8 > $out_vcf
bgzip -f $out_vcf



