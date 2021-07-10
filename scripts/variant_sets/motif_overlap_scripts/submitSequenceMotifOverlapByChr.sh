#!/bin/bash
##
##	submitMotifOverlapHocomocoByChr.sh
##
##	EDF 1/23/20
##



cd ~/data/gtex/v8/motif_overlap/

source ~/.bashrc



hocomoco_file=input_files/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv
line=`head -n $(($SLURM_ARRAY_TASK_ID + 1)) $hocomoco_file | tail -n 1`
TF=`echo $line | awk '{print $2}'`
qual=`echo $line | awk '{print $4}'`

final_file=motif_overlap_vars/$TF.HOCOMOCOv11_core_annotation_HUMAN_mono.$qual.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz
if [[ ! -e $final_file ]]
then
	num=22
	job_name=$TF\_MotifOverlapHocomocoByChr
	sbatch --job-name=$job_name --output=`pwd`/log/$job_name.%A.%a.out \
	        --partition=pe2 --array=1-$num --time=20:00:00 \
	        --mail-type=BEGIN,END,FAIL --mail-user=eflynn@nygenome.org\
	        scripts/getSequenceMotifOverlapByChr.sh $(($SLURM_ARRAY_TASK_ID + 1))
else
	echo $final_file already exists.
fi


