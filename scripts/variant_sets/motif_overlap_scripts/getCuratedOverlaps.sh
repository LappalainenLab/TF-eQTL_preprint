#!/bin/bash
##
##	getCuratedOverlaps.sh
##
##	- Select TFs that also have ChIPseq
##	- Select 5% MAF vars only
##	- Sum overlaps for those TFs
##
##	- Rewrote to work per TF
##
##	EDF 2/13/20
##

cd ~/data/gtex/v8/motif_overlap/
module load tabix

chr=$SLURM_ARRAY_TASK_ID

in_file=$1	## do for both change and either
in_header=`echo $in_file | sed 's/txt.gz/header/'`

per=$2		## maf percent cutoff
maf_gts=input_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF$per.noX.vcf.gz

out_file=split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt

echo in file is $in_file
echo maf is $per
echo maf file is $maf_gts
echo chr is $chr
echo out file is $out_file

#tfs_chipseq=(`awk '{if (NR>1) {print $1} }' input_files/TFs.info.txt`)
tfs_chipseq=input_files/TFs.info.txt

## select tfs that also have chipseq
keep_cols="1,2,3,4,5"
i=6
for col in `cat $in_header | cut -f 6-`
do
	if [[ "`grep -P ^$col"\t" $tfs_chipseq`" ]]
	then
		keep_cols+=","$i
	fi
	let i+=1
done


## cycle through maf 5% vars
#cat $in_header | cut -f $keep_cols | awk '{OFS="\t"; print $0,"sum"}' > $out_file
echo "Output columns will be..."
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
			for (i=6; i<=NF; i++) {sum+=$i} 
			print $0,sum }'
done >> $out_file


## prepare final output file
#source ~/.bashrc
#module load tabix
#bgzip -f $out_file
#tabix -s 2 -b 3 -e 3 $out_file.gz


echo "Done."
date

