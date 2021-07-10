#!/bin/bash
##
##	correlateTFs_wrapper.sh
##
##	Script to run correlateTFs.R
##
##	EDF 4/29/20
##

cd ~/projects/crosstiss_tf_corrs/

tf_file=input_files/TFs.info.curated.txt

maf=$1
tf=$2
chr=$SLURM_ARRAY_TASK_ID

out_file=correlations/by_tf/by_chr/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.$tf.chr$chr.txt
if [[ -e $out_file ]] && [[ `wc -l $out_file | awk '{print $1}'` -gt 11 ]]
then
	echo $out_file already exists and has more than eleven lines
	echo Exiting...
else
	echo $out_file does not exist or has eleven or fewer lines
	echo Running correlations...

	module load R/3.6.0
	Rscript scripts/correlateTFs_med0_bychr.R $maf $tf $chr
fi

echo "Done."
date



