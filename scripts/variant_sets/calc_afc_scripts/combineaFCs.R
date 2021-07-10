#!/usr/bin/Rscript
##
##  combineaFC_bychr.sh
##
##  Script to combine aFC calcs
##    across tissues
##
##  EDF 4/22/2020
##

library(dplyr)
library(tidyr)

setwd("~/projects/crosstiss_tf_corrs/afcs/")
args = commandArgs(trailingOnly=TRUE)
maf=args[1]
#maf="05"
chr=args[2]
#chr=20

tiss_trans = read.table("../input_files/tissue_translation_colors_v8.txt",
                        header=TRUE,sep='\t')

print("Combining aFCs with covariates...")
tiss_files = list.files(path=list.dirs("by_tiss"),
                        pattern=paste0("*.afcs.cov.overlap_vars.MAF",maf,".chr",chr,".txt"),
                        full.names=TRUE)
tiss_data = lapply(tiss_files, function(file) {
  tiss_long=strsplit(file,"/")[[1]][2]
  tiss_short=tiss_trans %>% filter(TISSUE_NAME==tiss_long) %>% pull(TISSUE_ABBRV) %>% as.character()
  print(c(tiss_long, tiss_short))
  read.table(file, header=TRUE, sep='\t') %>%
    select(sid,pid,!!tiss_short:="log2_aFC")
})

all_tiss_table = Reduce( function(x,y) {
    merge(x,y, by=c('sid','pid'), all=TRUE) },
    tiss_data) %>%
  separate(sid, c("chr","pos"), sep="_", remove=FALSE) %>%
  select(sid, chr, pos, pid,
         sort(everything()[-c(1,2,3,4)])) %>%
  arrange(as.numeric(pos))

write.table(all_tiss_table,
            paste0("combined/by_chr/all_tiss.afcs.cov.overlap_vars.MAF",maf,".chr",chr,".txt"),
            col.names = TRUE, row.names = FALSE,
            sep='\t', quote=FALSE)




## NOW FOR NO COVARIATES DATA
print("Combining aFCs without covariates...")

tiss_files = list.files(path=list.dirs("by_tiss"),
                        pattern=paste0("*.afcs.nocov.overlap_vars.MAF",maf,".chr",chr,".txt"),
                        full.names=TRUE)
tiss_data = lapply(tiss_files, function(file) {
  tiss_long=strsplit(file,"/")[[1]][2]
  tiss_short=tiss_trans %>% filter(TISSUE_NAME==tiss_long) %>% pull(TISSUE_ABBRV) %>% as.character()
  print(c(tiss_long, tiss_short))
  read.table(file, header=TRUE, sep='\t') %>%
    select(sid,pid,!!tiss_short:="log2_aFC")
})

all_tiss_table = Reduce( function(x,y) {
  merge(x,y, by=c('sid','pid'), all=TRUE) },
  tiss_data) %>%
  separate(sid, c("chr","pos"), sep="_", remove=FALSE) %>%
  select(sid, chr, pos, pid,
         sort(everything()[-c(1,2,3,4)])) %>%
  arrange(as.numeric(pos))

write.table(all_tiss_table,
            paste0("combined/by_chr/all_tiss.afcs.nocov.overlap_vars.MAF",maf,".chr",chr,".txt"),
            col.names = TRUE, row.names = FALSE,
            sep='\t', quote=FALSE)
