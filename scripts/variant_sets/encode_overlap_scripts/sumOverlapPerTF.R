#!/usr/bin/Rscript
##
##  calcPeakPerTF.R
##
##  Given a Peak overlap file, sum the number of 
##    overlaps per TF and DNase
##



setwd("~/data/gtex/v8/encode_overlap")
#args = commandArgs(trailingOnly=TRUE)
in_file = "ENCODE_TF_overlap.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt.gz"
out_file = "ENCODE_TF_overlap_by_TF.GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.txt"

## read in TFBS data
TFBS_header = readLines(in_file,n=1)

#TFBS_counts = read.table("ENCODE_GRCh38/TfDNase_split/TfDNase_split.aa")
#TFBS_header = read.table("ENCODE_GRCh38/WGS_Peak_overlap_TfDNase.header",
#                         header=TRUE,comment.char='')
#names(TFBS_counts) <- names(TFBS_header)


TFs = sort(unique( sapply( names(TFBS_counts)[-c(seq(1,3))] , 
                      function(str) {strsplit(str,'[_]')[[1]][1]} ) ))

TF_counts = as.data.frame( sapply(TFs, function(tf) {
  tf_cols = grep(paste('^',tf,'_',sep=''),names(TFBS_counts))
  apply(TFBS_counts[tf_cols], 1, function(row) {
    sum(row)
  })
}) )


out_TFBS_counts = cbind(TFBS_counts[c(seq(1,3))],
                        TF_counts)
write.table(out_TFBS_counts,file=out_file,
            col.names=TRUE,row.names=FALSE,
            quote=FALSE,sep='\t')




