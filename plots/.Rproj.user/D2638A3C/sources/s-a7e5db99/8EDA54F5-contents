#!/usr/bin/Rscript
##
##  fig5_APBB1IP.R
##
##  EDF 6/9/2021
##

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(ggsci)
library(scales)
library(ggrepel)

setwd("~/projects/MANUSCRIPT/")

gwas_eqtl_pvals = read.table("data_tables/GxE_GWAS/APBB1IP_gwas_eqtl_pvals.txt",
                             header=TRUE,sep='\t',
                             stringsAsFactors = FALSE)
overlap_var = "chr10_26438309_T_A_b38"

gwas_eqtl_pvals %>%
  filter(data_type %in% c('Astle_et_al_2016_Eosinophil_counts','Thyroid')) %>%
  mutate(data = ifelse(data_type == 'Astle_et_al_2016_Eosinophil_counts', 'Eosinophil counts GWAS',
                       ifelse(data_type == 'Thyroid', 'Thyroid APBB1IP eQTL', NA))) %>%
  mutate(r2_cat = factor(ifelse(variant_id == overlap_var, 'lead',
                                ifelse(r2>.8,'.8<r2<=1',
                                       ifelse(r2>.6,'.6<r2<=.8',
                                              ifelse(r2>.4,'.4<r2<=.6',
                                                     ifelse(r2>.2,'.2<r2<=.4','r2<=.2'))))),
                         levels=c('.8<r2<=1','.6<r2<=.8','.4<r2<=.6','.2<r2<=.4','r2<=.2','lead'))) %>%
  arrange(r2) %>%
  arrange(label) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r2_cat,
                 shape=r2_cat=='lead',
                 size=label),
             stroke=.1) +
  geom_point(aes(x=ifelse(label=='special',pos,NA),
                 fill=r2_cat,
                 shape=r2_cat=='lead',
                 size=label),
             stroke=.5) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(21,23)) +
  scale_size_manual(values=c(1,2)) +
  facet_wrap(~data,
             scales = 'free_y',
             nrow=3) +
  xlim(posi-10^5,posi+10^5)
ggsave("plots/fig5_gwas_eqtl_pvals.pdf",
       height=4,width=6)




thy_data = read.table("data_tables/GxE_GWAS/APBB1IP_thyroid_data.txt",
                      header=TRUE,sep='\t')
thy_corr = read.table("data_tables/GxE_GWAS/APBB1IP_thyroid_corr.txt",
                      header=TRUE,sep='\t')
thy_data %>%
  ggplot(aes(log10(gene_expr),log10(tf_expr))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=thy_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)'))
ggsave("plots/fig5_APBB1IP_tfieqtl_thy.pdf",
       height=3,width=4)



pit_data = read.table("data_tables/GxE_GWAS/APBB1IP_pituitary_data.txt",
                      header=TRUE,sep='\t')
pit_corr = read.table("data_tables/GxE_GWAS/APBB1IP_pituitary_corr.txt",
                      header=TRUE,sep='\t')
pit_data %>%
  ggplot(aes(log10(gene_expr),log10(tf_expr))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=pit_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)'))
ggsave("plots/fig5_APBB1IP_tfieqtl_pit.pdf",
       height=3,width=4)


