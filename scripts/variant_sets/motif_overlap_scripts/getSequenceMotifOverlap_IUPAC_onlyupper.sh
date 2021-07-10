#!/bin/bash
##
##	calcMotifOverlap_IUPAC_onlyupper.sh
##
##	EDF 6/15/2018
##

cd ~/data/gtex/v8/motif_overlap/

source ~/.bashrc
module load samtools
module load tabix

#TF=$1
#motif=$2
hocomoco_file=input_files/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv
line=`head -n $(($1 + 1)) $hocomoco_file | tail -n 1`
TF=`echo $line | awk '{print $2}'`
qual=`echo $line | awk '{print $4}'`
motif=`echo $line | awk '{print $6}'`
#motif=`echo $motif | sed 's/^[a-z]*//' | sed 's/[a-z]$//'`
## since we're not using lowercase ends, just remove them so it's shorter/more efficient
## nevermind, still interesting if var falls in region that was identified by motif even if N

#CaVEMaN_sp_med0_file='CaVEMaN/sp_med0.TFpeak_overlap.txt.gz'
input_vcf=input_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz
SNP_col=3

REF=~/data/GTEx_Analysis_v8/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
#out_file=HOCOMOCO_consensus/CaVEMaN_sp_med0.uniq_snps.$TF\.onlyuppermotif.txt
out_file=motif_overlap_vars/$TF.`basename $hocomoco_file | sed 's/tsv//'`$qual.onlyuppermotif.`basename $input_vcf | sed 's/vcf.gz//'`txt


echo "Testing TF $TF motif $motif, quality $qual, for overlap with $input_vcf"


regex_motif=`echo $motif | sed 's/[a-z]/N/g' | sed 's/R/[AG]/g' | sed 's/Y/[CT]/g' | \
	sed 's/S/[GC]/g' | sed 's/W/[AT]/g' | sed 's/K/[GT]/g' | sed 's/M/[AC]/g' | \
	sed 's/B/[CGT]/g' | sed 's/D/[AGT]/g' | sed 's/H/[ACT]/g' | sed 's/V/[ACG]/g' | \
	sed 's/N/[ACGT]/g'`
b_regex_motif=`echo $motif | sed 's/[a-z]/N/g' | rev | sed 's/R/[AG]/g' | sed 's/Y/[CT]/g' | \
	sed 's/S/[GC]/g' | sed 's/W/[AT]/g' | sed 's/K/[GT]/g' | sed 's/M/[AC]/g' | \
	sed 's/B/[CGT]/g' | sed 's/D/[AGT]/g' | sed 's/H/[ACT]/g' | sed 's/V/[ACG]/g' | \
	sed 's/N/[ACGT]/g' | sed -e 's/A/t/g' -e 's/T/a/g' -e 's/C/g/g' -e 's/G/c/g' | \
	awk '{print toupper($0)}' `
let mot_len=`echo $motif | wc -c`-1		## wc -c adds 1
echo -e "Motif: $motif \nRegex motif: $regex_motif; Backwards motif: $b_regex_motif"


echo snp chr pos ref alt ref_mot alt_mot either change | tr ' ' '\t' > $out_file

gunzip < $input_vcf | grep -v "#" | while read line
do
	snp=`echo $line | cut -d ' ' -f $SNP_col`
	#echo $SNP
	chr=`echo $snp | cut -d_ -f 1`
	pos=`echo $snp | cut -d_ -f 2`
	ref=`echo $snp | cut -d_ -f 3`
	alt=`echo $snp | cut -d_ -f 4`

	let sta=$pos-$mot_len+1
	let end=$pos+`echo $ref | wc -c`-1+$mot_len-2

	#echo $pos $ref $alt $sta $end

	seq=`samtools faidx $REF $chr':'$sta'-'$end | tail -n +2`
	mid1=$mot_len-1		
	## 0-based strings
	let mid2=$mot_len-1+`echo $ref | wc -c`-1
	alt_seq=`echo ${seq:0:$mid1}``echo $alt``echo ${seq:$mid2}`

	#echo $seq $alt_seq

	if [[ "`echo $seq | grep -e $regex_motif -e $b_regex_motif`" ]]; then ref_mot=Y; else ref_mot=N; fi
	if [[ "`echo $alt_seq | grep -e $regex_motif -e $b_regex_motif`" ]]; then alt_mot=Y; else alt_mot=N; fi
	if [[ $ref_mot == 'Y' || $alt_mot == 'Y' ]]; then either=Y; \
		if [[ $ref_mot == 'N' || $alt_mot == 'N' ]]; then change=Y; else change=N; fi; \
		else either=N; change=N; fi

	#echo $ref_mot $alt_mot

	echo $snp $chr $pos $ref $alt $ref_mot $alt_mot $either $change | tr ' ' '\t'

done >> $out_file

bgzip $out_file
tabix -s 2 -b 3 -e 3 $out_file.gz


echo 'Done'
date








