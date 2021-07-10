#!/bin/bash


cd ~/projects/crosstiss_tf_corrs/

out_file=correlations/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.all.top.adjp2.txt

echo tf gene var sp_rho sp_p n_var meff sp_p_meff sp_p_meff_bh_tf | \
	tr ' ' '\t' > $out_file
for file in correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.*.top.adjp.txt
do
	ls $file
	tf=`basename $file | awk -F[.] '{print $5}'`

	cat $file | tail -n +2 | \
		awk -v var=$tf '{OFS="\t"; print var,$1,$7,$8,$4,$2,$3,$5,$12}' \
		>> $out_file
done

gzip $out_file


