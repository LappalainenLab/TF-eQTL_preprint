#!/bin/bash
##
##	sumOverlapByTF.sh
##
##	EDF 2/3/2020
##

cd ~/data/gtex/v8/encode_overlap
orig_file="ENCODE_TF_overlap.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.txt.gz"
temp_num=`echo $(($SLURM_ARRAY_TASK_ID-1)) | sed 's/^/00/'`
in_file=overlap_split/split.${temp_num: -3}
#out_file="ENCODE_TF_overlap_by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt"
out_file=`echo $in_file | sed 's/$/.out/'`

echo "In file is $in_file"
echo "Out file will be $out_file"

TFBS_header=`zcat $orig_file | head -n 1`
#all_cols=`echo $TFBS_header | tr ' ' '\n'`
TFs=(`echo $TFBS_header | tr ' ' '\n' | tail -n +4 | awk -F_ '{print $1}' | sort | uniq`)
TF_cols=(`for tf in ${TFs[@]}
	do 
		echo $TFBS_header | tr ' ' '\n' | grep -n $tf | 
			awk -F':' '{print $1}' | tr '\n' '_' | sed 's/_$/\n/'
	done`)


#echo `echo $TFBS_header | tr ' ' '\n' | head -n 3` ${TFs[@]} all | \
#	tr ' ' '\t' > $out_file
while read line
do
	#echo $line
	out=`echo $line | awk '{print $1,$2,$3}'`
	for cols in ${TF_cols[@]}
	do
		#echo $cols
		out+=" "$(echo $line | cut -d ' ' -f `echo $cols | tr '_' ','` | \
			tr ' ' '+' | bc) 
	done
	echo $out | awk '{OFS="\t"; count=0;
			for (i=4; i<=NF; i++) {if ($i>0) {$i=1; count+=1} }
			print $0,count }'
done < $in_file > $out_file

