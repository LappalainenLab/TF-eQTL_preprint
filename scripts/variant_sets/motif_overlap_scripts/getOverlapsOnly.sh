#!/bin/bash
##
## getOverlapsOnly.sh
##
## EDF 4/8/2020
##

set=$1
per=$2
in_file=overlaps/HOCOMOCO_TF_motif_overlap.$set.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.curated_set.MAF$per.txt.gz
out_file=overlaps/HOCOMOCO_TF_motif_overlap.$set.GTEx_v8.curated_set.MAF$per.txt
out_vars=overlaps/HOCOMOCO_TF_motif_overlap.$set.GTEx_v8.curated_set.MAF$per.vars.list
zcat $in_file | awk '{if ($NF > 0) {print $0} }' | \
        tee $out_file | \
        awk '{if (NR > 1) {print $1} }' > $out_vars
module load tabix
bgzip -f $out_file
tabix -s 2 -b 3 -e 3 $out_file.gz

