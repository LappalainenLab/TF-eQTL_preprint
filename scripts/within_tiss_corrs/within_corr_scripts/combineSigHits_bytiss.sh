#!/bin/bash
##
##	combineSigHits_bytf.sh
##
##	EDF 7/9/2020
##


tiss=$1

echo $tiss
num_files=`ls tensorqtl/$tiss/*/*.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | wc -l | awk '{print $1}'`

all_out_file=tensorqtl/$tiss.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz
out_file1=tensorqtl/$tiss.norm.ieqtl.sig_assoc.fdr05.txt
out_file2=tensorqtl/$tiss.norm.ieqtl.sig_assoc.fdr20.txt

zcat tensorqtl/Adipose_Subcutaneous/AHR/Adipose_Subcutaneous.AHR.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | \
        awk '{OFS="\t"; if (NR == 1) {print $0,"tf","tiss"} }' > $all_out_file
zcat tensorqtl/Adipose_Subcutaneous/AHR/Adipose_Subcutaneous.AHR.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | \
	awk '{OFS="\t"; if (NR == 1) {print $0,"tf","tiss"} }' > $out_file1
zcat tensorqtl/Adipose_Subcutaneous/AHR/Adipose_Subcutaneous.AHR.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz | \
	awk '{OFS="\t"; if (NR == 1) {print $0,"tf","tiss"} }' > $out_file2

#for tiss in `cat input_files/representative_tissues.20.txt | awk '{ if (NR > 1) {print $1} }'`
for file in tensorqtl/$tiss/*/*.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz
do
	tf=`echo $file | awk -F'/' '{print $3}'`
	echo $tf
	zcat $file | \
		awk -v var1=$tf -v var2=$tiss '{OFS="\t"; print $0, var1, var2 }' | \
		tee -a $all_out_file | \
		awk '{OFS="\t"; if ($18 < 0.20) {print $0} }' | \
		tee -a $out_file2 | \
		awk '{OFS="\t"; if ($18 < 0.05) {print $0} }' | \
		tee -a $out_file1 | \
		wc -l
done

echo "Done."
date
