#!/bin/bash
##
##	combineCuratedOverlaps.sh
##
##
##	EDF 2/13/20
##

cd ~/data/gtex/v8/motif_overlap/
module load tabix

in_file=$1	## do for both change and either
in_header=`echo $in_file | sed 's/txt.gz/header/'`
#out_file=split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt

per=$2		## maf percent cutoff
maf_gts=input_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF$per.noX.vcf.gz

out_file=`basename $in_file | sed 's/txt.gz//'`curated_set.MAF$per.txt

echo orig file is $in_file
echo maf is $per
echo maf file was $maf_gts
echo out file is $out_file

## print new header
#tfs_chipseq=(`awk '{if (NR>1) {print $1} }' input_files/TFs.info.txt`)
tfs_chipseq=input_files/TFs.info.txt
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
cat $in_header | cut -f $keep_cols | awk '{OFS="\t"; print "#"$0,"sum"}' > $out_file

## cycle through chr
for chr in `seq 1 22`
do
	this_file=split/`basename $in_file | sed 's/txt.gz//'`chr$chr.curated_set.MAF$per.txt
	cat $this_file
done >> $out_file


## prepare final output file
#source ~/.bashrc
#module load tabix
wc -l $out_file
bgzip -f $out_file
tabix -s 2 -b 3 -e 3 $out_file.gz


echo "Done."
date

