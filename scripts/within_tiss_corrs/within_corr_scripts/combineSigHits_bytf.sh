#!/bin/bash
##
##	combineSigHits_bytf.sh
##
##	EDF 7/9/2020
##


tf=$1

echo $tf
num_files=`ls tensorqtl/*/$tf/*.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | wc -l | awk '{print $1}'`

if [[ $num_files -eq 200 ]];
then
	out_file1=tensorqtl/$tf.norm.ieqtl.sig_assoc.fdr05.txt
	out_file2=tensorqtl/$tf.norm.ieqtl.sig_assoc.fdr20.txt
else
	out_file1=tensorqtl/partial_tf_files/$tf.norm.ieqtl.sig_assoc.fdr05.partial.txt
	out_file2=tensorqtl/partial_tf_files/$tf.norm.ieqtl.sig_assoc.fdr20.partial.txt
fi

zcat tensorqtl/Adipose_Subcutaneous/AHR/Adipose_Subcutaneous.AHR.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | \
	awk '{OFS="\t"; if (NR == 1) {print $0,"tf","tiss"} }' > $out_file1
zcat tensorqtl/Adipose_Subcutaneous/AHR/Adipose_Subcutaneous.AHR.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | \
	awk '{OFS="\t"; if (NR == 1) {print $0,"tf","tiss"} }' > $out_file2

#for tiss in `cat input_files/representative_tissues.20.txt | awk '{ if (NR > 1) {print $1} }'`
for file in tensorqtl/*/$tf/*.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz
do
	tiss=`echo $file | awk -F'/' '{print $2}'`
	echo $tiss
	zcat $file | \
		awk -v var1=$tf -v var2=$tiss '{OFS="\t"; if ($18 < 0.20) {print $0, var1, var2} }' | \
		tee -a $out_file2 | \
		awk -v var1=$tf -v var2=$tiss '{OFS="\t"; if ($18 < 0.05) {print $0} }' | \
		tee -a $out_file1 | \
		wc -l
done

echo "Done."
date
