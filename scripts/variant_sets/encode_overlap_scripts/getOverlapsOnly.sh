#!/bin/bash
##
## getOverlapsOnly.sh
##
## EDF 4/9/2020
##

per=$1
in_file=overlaps/ENCODE_TF_overlap.by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.curated_set.MAF$per.txt.gz
out_file=overlaps/ENCODE_TF_ChIPseq_overlap.GTEx_v8.curated_set.MAF$per.txt
out_vars=overlaps/ENCODE_TF_ChIPseq_overlap.GTEx_v8.curated_set.MAF$per.vars.list
zcat $in_file | awk '{if ($NF > 0) {print $0} }' | \
        tee $out_file | \
        awk '{if (NR > 1) {print $3} }' > $out_vars
module load tabix
bgzip -f $out_file
tabix -s 1 -b 2 -e 2 $out_file.gz

