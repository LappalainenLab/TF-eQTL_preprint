#!/bin/bash
##
##	getCuratedOverlaps.sh
##
##	- Select TFs that also have ChIPseq
##	- Select 5% MAF vars only
##	- Sum overlaps for those TFs
##
##
##	EDF 2/13/20
##

cd ~/data/gtex/v8/encode_overlap/
module load tabix

chr=$SLURM_ARRAY_TASK_ID
per=$1

in_file=overlaps/ENCODE_TF_overlap.by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.txt.gz
in_header=`echo $in_file | sed 's/txt.gz/header/'`

maf_gts=genetic_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF$per.noX.vcf.gz

out_file=overlap_split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt

echo in file is $in_file
echo chr is $chr
echo maf is $per
echo maf file is $maf_gts
echo out file will be $out_file

#tfs_chipseq=(`awk '{if (NR>1) {print $1} }' input_files/TFs.info.txt`)
tfs_motif=../motif_overlap/input_files/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv

## select tfs that also have motif
keep_cols="1,2,3"
i=4
for col in `cat $in_header | cut -f 4-`
do
	if [[ "`grep -P "\t"$col"\t" $tfs_motif`" ]]
	then
		keep_cols+=","$i
	fi
	let i+=1
done


## cycle through maf 5% vars
#cat $in_header | cut -f $keep_cols | awk '{OFS="\t"; print $0,"sum"}' > $out_file
echo "output columns will be..."
cat $in_header | cut -f $keep_cols | awk '{OFS="\t"; print $0,"sum"}'
gunzip < $maf_gts | grep -P ^chr$chr"\t" | while read line
do
	chr=`echo $line | awk '{print $1}'`
	pos=`echo $line | awk '{print $2}'`
	snp=`echo $line | awk '{print $3}'`
	#echo $snp
	
	## sum new tfs
	tabix $in_file $chr:$pos-$pos | grep $snp | \
		cut -f $keep_cols | \
		awk '{OFS="\t"; sum=0; 
			for (i=4; i<=NF; i++) {sum+=$i} 
			print $0,sum }' 
done > $out_file


## prepare final output file
source ~/.bashrc
module load tabix
bgzip -f $out_file
tabix -s 1 -b 2 -e 2 $out_file.gz


echo "Done."
date

