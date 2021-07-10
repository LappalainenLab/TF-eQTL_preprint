#!/bin/bash
##
##	combineTopTFCorrs.sh
##
##	EDF 6/17/20
##

cd ~/projects/crosstiss_tf_corrs/

out_file=correlations/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.all.top.adjp.txt


echo tf `head -n 1 correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.ATF3.top.adjp.txt` | tr ' ' '\t' > $out_file

for file in correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.*.top.adjp.txt
do
	wc -l $file
	tf=`echo $file | awk -F. '{print $5}'`
	awk -v var=$tf '{OFS="\t";
		if (NR>1) {print var,$0} }' $file >> $out_file
done
bgzip $out_file

echo "Done."
date


