#!/usr/bin/Rscript
##
##  getTopProtVars.R
##
##  EDF 6/22/21
##

setwd("~/projects/overlap/")

library(dplyr)
library(tidyr)

prot_hits = read.table("replication/sig_assoc.fdr20.protein_hits_filtered.txt",
                       header=TRUE,sep='\t')
all_prot_corrs = read.table("../prot_corrs/sig_prot_corrs.fdr20.txt",
                            header=TRUE,sep='\t')

prot_hits_var = merge(prot_hits,
      all_prot_corrs %>% select(tf,gene,var.min),
      by.x=c('tf','phenotype_id'),
      by.y=c('tf','gene'))
write.table(prot_hits_var,
            "replication/sig_assoc.fdr20.protein_hits_filtered.vars.txt",
            quote=FALSE,sep='\t',
            col.names=TRUE,row.names=FALSE)


