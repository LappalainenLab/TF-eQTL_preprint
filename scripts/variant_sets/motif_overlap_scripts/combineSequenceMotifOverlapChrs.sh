#!/bin/bash
##
##	combineSequenceMotifOverlapHocomocoChrs.sh
##
##	EDF 1/23/20
##



cd ~/data/gtex/v8/motif_overlap/

source ~/.bashrc

hocomoco_file=input_files/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv
line=`head -n $(($SLURM_ARRAY_TASK_ID + 1)) $hocomoco_file | tail -n 1`
TF=`echo $line | awk '{print $2}'`
qual=`echo $line | awk '{print $4}'`

out_file=motif_overlap_vars/$TF.HOCOMOCOv11_core_annotation_HUMAN_mono.$qual.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt

if [[ "`ls motif_overlap_vars/$TF.chr*.HOCOMOCOv11_core_annotation_HUMAN_mono.$qual.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz`" ]]
then

	echo "Combining all chrs for $TF to $out_file"

	zcat motif_overlap_vars/$TF.chr10.HOCOMOCOv11_core_annotation_HUMAN_mono.$qual.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz | \
		head -n 1 | tr ' ' '\t' > $out_file

	for chr in `seq 1 22`
	do
		zcat motif_overlap_vars/$TF.chr$chr.HOCOMOCOv11_core_annotation_HUMAN_mono.$qual.onlyuppermotif.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz | \
			tail -n +2 >> $out_file
	done

	known_length=10403249
	if [[ "`wc -l $out_file | awk '{print $1}'`" -ne 10403249 ]]
	then
		echo "ERROR!! FILE IS NOT CORRECT LENGTH"
		wc -l $out_file
		mv motif_overlap_vars/$TF.chr*txt motif_overlap_vars/by_chr_bad/
		mv motif_overlap_vars/$TF.chr*txt.gz motif_overlap_vars/by_chr_bad/
	else
		echo bgzipping and tabixing file
		bgzip -f $out_file
		tabix -s 2 -b 3 -e 3 -S 1 $out_file.gz
		mv motif_overlap_vars/$TF.chr*txt.gz motif_overlap_vars/by_chr/	
	fi
else
	echo "No input files for $TF"
	echo "Does output file exist?"
	ls -l $out_file.gz
fi
	
echo "Done."
date

