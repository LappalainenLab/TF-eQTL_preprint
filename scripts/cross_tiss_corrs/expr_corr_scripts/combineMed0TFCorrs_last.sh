#!/bin/bash
##
## combineTFCorrs_05.sh
##
##	EDF 5/5/20
##


cd ~/projects/crosstiss_tf_corrs/

#tf_file=input_files/TFs.info.curated.txt

maf=$1
#tf=`head -n $(echo $SLURM_ARRAY_TASK_ID"+1" | bc) $tf_file | \
#        tail -n 1 | awk '{print $1}'`

#echo "TF is $tf"


## combine all tiss files
echo "Combining all tiss files..."

out_file=correlations/cross_tiss_tf_expr_corrs.med0.curated_set.MAF$maf.txt

zcat correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF$maf.ATF3.txt.gz | \
	head -n 1 > $out_file
for file in correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF$maf.*.txt.gz
do
	tf=`basename $file | awk -F'.' '{print $5}'`
	echo -n $tf" "
	zcat $file | \
		tail -n +2 | \
		tee -a $out_file | \
		wc -l
done

wc -l $out_file
bgzip -f $out_file
## not in order by var, so can't tabix....



echo "done."
date

