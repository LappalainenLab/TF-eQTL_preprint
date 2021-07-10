#!/bin/bash
##
##	combineProtSigHits_fdr20.sh
##
##	EDF 9/17/2020
##

cd ~/projects/prot_corrs/

out_file=sig_prot_corrs.fdr05.txt
head -n 1 crosstiss_tf_corrs/correlations/by_tf/stats/cross_tiss_tf_protein_corrs.med0.MAF05.AHR.top.adjp.txt | \
	awk '{OFS="\t"; print $0,"tf"}' > $out_file

for file in crosstiss_tf_corrs/correlations/by_tf/stats/cross_tiss_tf_protein_corrs.med0.MAF05.*.top.adjp.txt
do
	echo $file
	tf=`basename $file | awk -F'.' '{print $4}'`
	awk -v var=$tf '{ if ($10 < 0.05) {OFS="\t"; print $0,var} }' $file >> $out_file
done

echo Done.
date



