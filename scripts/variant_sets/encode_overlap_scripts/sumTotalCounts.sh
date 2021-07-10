#!/bin/bash
##
##	sumTotalCounts.sh
##
##	EDF 2/11/20
##

cd ~/data/gtex/v8/encode_overlap/

temp_num=`echo ${SLURM_ARRAY_TASK_ID-1} | sed 's/^/00/'`
in_file=overlap_split/split.${temp_num: -3}.out
#out_file=ENCODE_TF_overlap.by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt
out_file=`echo $in_file | sed 's/$/.2/'`

echo "In file is $in_file"
echo "Out file will be $out_file"

#zcat ENCODE_TF_overlap.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz | \
#	head -n 1 | awk '{OFS="\t"; print $0,"count"}' > $out_file

#files=`ls overlap_split/split.*.out`
#cat $files | \
#	awk '{OFS="\t"; count=0; for (i=4; i<=NF; i++) {if ($i>0) {count+=1} } print $0,count}' \
#	>> $out_file

cat $in_file | awk '{OFS="\t"; count=0;
	for (i=4; i<=NF; i++) { if ($i>0) { count+=1; $i=1 } }
	print $0,count}' \
	> $out_file



module load tabix
bgzip $out_file
tabix -s 1 -b 2 -e 2 -S 1 $out_file.gz

echo "Done."
date


