#!/bin/bash
##
##	combineCuratedOverlaps.sh
##
##
##	EDF 2/13/20
##

cd ~/data/gtex/v8/encode_overlap/
module load tabix

per=$1

in_file=overlaps/ENCODE_TF_overlap.by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF01_GTonly.txt.gz
in_header=`echo $in_file | sed 's/txt.gz/header/'`
#out_file=split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt

maf_gts=genetic_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF$per.noX.vcf.gz

out_file=overlaps/`basename $in_file | sed 's/txt.gz//'`curated_set.MAF$per.txt

echo orig file is $in_file
echo maf is $per
echo maf file was $maf_gts
echo out file is $out_file

## print new header
#tfs_chipseq=(`awk '{if (NR>1) {print $1} }' input_files/TFs.info.txt`)
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
cat $in_header | cut -f $keep_cols | awk '{OFS="\t"; print "#"$0,"sum"}' > $out_file

## cycle through chr
for chr in `seq 1 22`
do
	this_file=overlap_split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt
	zcat $this_file
done >> $out_file


## prepare final output file
#source ~/.bashrc
module load tabix
wc -l $out_file
bgzip -f $out_file
tabix -s 1 -b 2 -e 2 $out_file.gz


echo "Done."
date

