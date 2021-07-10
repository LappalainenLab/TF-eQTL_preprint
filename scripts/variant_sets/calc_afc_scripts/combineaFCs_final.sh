#!/bin/bash
##
##	combineaFCs_final.sh
##
##	Last step for combining aFCs
##
##	EDF 4/23/20
##

maf=$1
cov=$2
out_file=combined/all_tiss.afcs.$cov.overlap_vars.MAF$maf.txt

echo out file will be $out_file

head -n 1 combined/by_chr/all_tiss.afcs.$cov.overlap_vars.MAF$maf.chr19.txt > $out_file
for chr in `seq 1 22`
do
	tail -n +2 combined/by_chr/all_tiss.afcs.$cov.overlap_vars.MAF$maf.chr$chr.txt >> $out_file
done

wc -l $out_file
bgzip -f $out_file
tabix -s 2 -b 3 -e 3 -S 1 $out_file.gz

echo Done.
date








