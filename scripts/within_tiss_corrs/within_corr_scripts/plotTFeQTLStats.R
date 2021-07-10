#!/usr/bin/Rscript
##
##  plotTFeQTLStats.R
##
##  EDF 7/2/11
##

library(dplyr)
library(ggplot2)

setwd("~/projects/TFi-eQTL/")

tiss_info = read.table('input_files/tissue_info.txt',
                       header=TRUE,sep='\t')

all_tfeqtl = read.table("tensorqtl/all_tiss_ctalso.norm.ieqtl.sig_assoc.fdr20.txt",
                        header=TRUE,sep='\t')
sig_tfeqtl = read.table("eqtl_sig/all_tiss_ctalso.fdr20.sig_eqtls.txt",
                        header=TRUE,sep='\t') %>%
  unite(tiss_gene, tiss, gene, remove=FALSE) %>%
  unique()

table(all_tfeqtl$tiss)
table(sig_tfeqtl$tiss)

all_tfeqtl %>%
  unite(tiss_gene, tiss, phenotype_id, remove=FALSE) %>%
  mutate(is_sig = tiss_gene %in% as.character(sig_tfeqtl$tiss_gene)) %>%
  ggplot(aes(tiss)) +
  geom_bar(aes(fill=is_sig)) +
  theme_classic() +
  scale_fill_manual(values=c('gray','blue')) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab('20% FDR TF-eQTLs') +
  xlab('Tissue')

sig_tfeqtl %>% unite(tisi, tiss_gene, tf) %>% pull(tisi) %>% as.character() %>% unique() %>% length()
sig_tfeqtl %>% unite(tisi, tf, gene) %>% pull(tisi) %>% as.character() %>% unique() %>% length()

temp = sig_tfeqtl %>%
  group_by(tiss) %>%
  summarize(n=n()) %>%
  merge(tiss_info, by.x='tiss', by.y='TISSUE_NAME')

temp %>%
  ggplot(aes(SAMPLE_COUNT,n)) +
  geom_point(color=as.character(temp$TISSUE_RCOL)) +
  geom_text(aes(label=ifelse(n>10000,as.character(SMTSD),''))) +
  theme_classic() +
  ylab('20% FDR TF-ieQTLs') +
  xlab('Tissue Sample Count')




sig_tfeqtl_final = read.table('eqtl_sig/all.fdr20.sig_eqtls.txt',
                              header=TRUE, sep='\t') %>%
  unite(tf_gene, tf, gene, remove=FALSE)

hist(table(as.character(sig_tfeqtl_final$tf_gene)),
     breaks=c(0,1,2,3,4,5,6,7,8,9,10,11),
     main='Histogram of tissues per TF-eQTL',
     xlab='Tissues',
     ylab='TFeQTLs')



table(as.character(sig_tfeqtl_final$tf_gene)) %>% length()
sum(table(as.character(sig_tfeqtl_final$tf_gene)) > 1)
