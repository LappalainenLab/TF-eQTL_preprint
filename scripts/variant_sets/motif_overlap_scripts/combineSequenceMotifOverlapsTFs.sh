#!/bin/bash
##
##	combineSequenceMotifOverlapsTFs.sh
##
##	2/3/2020 EDF
##

cd ~/lap_lab_folder/data/gtex/v8/motif_overlap/

out_file1=HOCOMOCO_TF_motif_overlap.either.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.txt
#out_file1=HOCOMOCO_TF_motif_overlap.change.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.txt

echo out file will be $out_file1
ex_file=motif_overlap_vars/ZSCAN31.HOCOMOCOv11_core_annotation_HUMAN_mono.C.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz
files=`ls motif_overlap_vars/*.HOCOMOCOv11_core_annotation_HUMAN_mono.*.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz`
tfs=`ls $files | awk -F/ '{print $2}' | awk -F. '{print $1}'`

echo `zcat $ex_file | head -n 1 | cut -f1,2,3,4,5` $tfs sum | tr ' ' '\t' > $out_file1
command="paste <(zcat $ex_file | cut -f1,2,3,4,5)"
for file in $files
do
	command=$command" <(zcat $file | cut -f8)"
	#command=$command" <(zcat $file | cut -f9)"
done
echo $command
eval $command | tail -n +2 | \
	tr 'N' '0' | tr 'Y' '1' | \
	awk '{OFS="\t"; sum=0; for (i=6; i<=NF; i++) {sum+=$i} print $0,sum}' \
	>> $out_file1

module load tabix
bgzip $out_file1
tabix -s 2 -b 3 -e 3 -S 1 $out_file1.gz

echo "Done."
date
