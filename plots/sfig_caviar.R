#!/usr/bin/Rscript
##
##  sfig_caviar.R
##
##
##

library(dplyr)
library(tidyr)
library(ggplot2)

cav_hits = read.table("~/projects/TFi-eQTL/variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t')

table(cav_hits$gene) %>%
  hist(main="Variants in fine-mapped set per gene", 
       xlab="Variants",
       ylab="Genes")
median(table(cav_hits$gene))

table(cav_hits$var) %>%
  hist(main="Genes associated with each variant", 
       xlab="Genes",
       ylab="Variants")
median(table(cav_hits$var))

cross_hits = read.table("~/projects/MANUSCRIPT/data_tables/cross_tiss_corrs/cross_corrs_fdr05.txt",
                        header=TRUE,sep='\t')
prot_hits = read.table("~/projects/overlap/replication/sig_assoc.fdr20.protein_hits_filtered.txt",
                       header=TRUE,sep='\t')
with_hits = read.table("~/projects/MANUSCRIPT/data_tables/within_tiss_corrs/win_corrs_fdr20.sig_eqtl.txt",
                       header=TRUE,sep='\t')
moc_hits = read.table("~/projects/MANUSCRIPT/data_tables/combine_corrs/multi_or_with_cross_corrs.txt",
                      header=TRUE,sep='\t')

cav_hits %>%
  group_by(gene) %>%
  summarize(n_var=n()) %>%
  mutate(cross = gene %in% cross_hits$gene,
         prot = gene %in% prot_hits$phenotype_id,
         with = gene %in% with_hits$gene,
         moc = gene %in% moc_hits$gene) %>%
  pivot_longer(cols=c(cross,prot,with,moc),
               values_to='corr') %>%
  mutate(name=factor(name,
                     levels = c('cross','prot','with','moc'))) %>%
  ggplot(aes(name,n_var)) +
  geom_boxplot(aes(col=corr)) +
  theme_classic() +
  xlab("Corr type") +
  ylab("Number of vars per gene")
  




