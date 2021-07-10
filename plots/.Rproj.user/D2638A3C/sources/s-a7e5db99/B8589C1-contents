#!/usr/bin/Rscript
##
##  fig2_plotExpl.R
##
##  EDF 5/4/2021
##

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(seqminer)


setwd("~/projects/MANUSCRIPT/")

tiss_info = read.table("../TFi-eQTL/input_files/tissue_info.txt",
                       header=TRUE, sep='\t') %>%
  mutate(tiss_type = factor(
    c(rep('fat',2),'organ',rep('epithelial',4),
      rep('brain',13),
      'fat','epithelial','blood',rep('epithelial',8),
      rep('muscle',2),rep('organ',3),
      'epithelial','organ','muscle','nervous',rep('organ',2),
      'nervous','organ',rep('epithelial',3),'blood',
      'epithelial','organ','nervous','epithelial','epithelial','blood')
  ),
  tiss_type_2 = factor(ifelse(tiss_type=='epithelial','epithelial',
                              ifelse(tiss_type=='nervous','nervous',
                                     'internal'))))




descriptioni="MS4A14"
tfi='FOSL2'
overlap_var='chr11_60389634_T_G_b38'

gene_info = read.table("../crosstiss_tf_corrs/input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                       header=TRUE,sep='\t')
genei = gene_info %>%
  filter(description==descriptioni) %>%
  pull(gene) %>%
  as.character()

cross_tiss_stats = read.table(paste0("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",
                                     tfi,".top.adjp.txt"),
                              header=TRUE,sep='\t') %>%
  filter(gene == genei)


sig_tiss_stats = read.table("../overlap/replication/tiss_filt/all.fdr20.sig_eqtls.txt",
                            header=TRUE,sep='\t') %>%
  filter(tf==tfi, gene == genei)
sig_tiss = sig_tiss_stats %>% pull(tiss) %>% as.character()
sig_tiss_short = tiss_info %>%
  filter(TISSUE_NAME == sig_tiss[1]) %>%
  pull(TISSUE_ABBRV)
ns_tiss_short = 'ADPVSC'


ct_tf_expr = read.table("../crosstiss_tf_corrs/input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                          header=TRUE,sep='\t') %>%
  filter(description == tfi) %>%
  select(-c(gene,description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='tiss_short')

ct_colnames = read.table("../crosstiss_tf_corrs/afcs/combined/all_tiss.afcs.cov.overlap_vars.MAF05.header")
ct_afc = tabix.read.table("../crosstiss_tf_corrs/afcs/combined/all_tiss.afcs.cov.overlap_vars.MAF05.txt.gz",
                          paste0(strsplit(overlap_var,"_")[[1]][1],":",
                                 strsplit(overlap_var,"_")[[1]][2],"-",
                                 strsplit(overlap_var,"_")[[1]][2])) %>%
  filter(V1 == overlap_var,
         V4 == genei) %>%
  select(-c(V1,V2,V3,V4)) %>%
  t() %>% as.data.frame() %>%
  mutate(tiss_short = ct_colnames[-c(1,2,3,4)] %>% t() %>% unlist())


ct_data = merge(ct_tf_expr, ct_afc,
                by='tiss_short') %>%
  merge(tiss_info, 
        by.x='tiss_short',
        by.y='TISSUE_ABBRV')

ct_cor = cor.test(ct_data$V1.x, ct_data$V1.y,
                  method='spearman')

ct_data %>%
  ggplot(aes(log10(V1.x),V1.y)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=tiss_type)) +
  theme_classic() +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab('eQTL aFC') +
  ggtitle(paste(descriptioni, genei, '\n', overlap_var)) +
  geom_label(x=.5,y=.1,
             label=paste0('rho=',round(ct_cor$estimate,2),'\n',
                          'p=',round(ct_cor$p.value,2)),
             hjust=0) +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('hotpink2','goldenrod2','deepskyblue3',
                              'chocolate2','slateblue3','olivedrab3',
                              'antiquewhite3'))

ct_data %>%
  ggplot(aes(log10(V1.x),V1.y)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=ifelse(tiss_short %in% c(as.character(sig_tiss_short),ns_tiss_short),
                            'zzz',tiss_type),
                 fill=tiss_type,
                 shape=(tiss_short %in% c(as.character(sig_tiss_short),ns_tiss_short)))) +
  theme_classic() +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab('eQTL aFC') +
  ggtitle(paste(descriptioni, genei, '\n', overlap_var)) +
  geom_label(x=.5,y=.1,
             label=paste0('rho=',round(ct_cor$estimate,2),'\n',
                          'p=',round(ct_cor$p.value,2)),
             hjust=0) +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('hotpink2','goldenrod2','deepskyblue3',
                             'chocolate2','slateblue3','olivedrab3',
                             'antiquewhite3')) +
  scale_color_manual(values=c('hotpink2','goldenrod2','deepskyblue3',
                              'chocolate2','slateblue3','olivedrab3',
                              'antiquewhite3','black')) +
  scale_shape_manual(values=c(16,21))
ct_data %>%
  ggplot(aes(log10(V1.x),V1.y)) +
  geom_point(aes(col=ifelse(tiss_short %in% c(as.character(sig_tiss_short),ns_tiss_short),
                            'zzz',tiss_type),
                 fill=tiss_type,
                 shape=(tiss_short %in% c(as.character(sig_tiss_short),ns_tiss_short))),
             size=2) +
  theme_classic() +
  xlab(paste0('log10(median ',tfi,' TPM)')) +
  ylab(paste0(descriptioni,' eQTL aFC')) +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('hotpink2','goldenrod2','deepskyblue3',
                             'chocolate2','slateblue3','olivedrab3',
                             'antiquewhite3')) +
  scale_color_manual(values=c('hotpink2','goldenrod2','deepskyblue3',
                              'chocolate2','slateblue3','olivedrab3',
                              'antiquewhite3','black')) +
  scale_shape_manual(values=c(16,21))
ggsave(paste(paste0('plots/',tfi),descriptioni,overlap_var,"cross.pdf",sep='.'),
       height=3, width=3)


tiss_tf_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short,".tpm.gct"),
                            header=TRUE, sep='\t',
                            skip = 2) %>%
  filter(Description == tfi) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

tiss_gene_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short,".tpm.gct"),
                            header=TRUE, sep='\t',
                            skip = 2) %>%
  filter(Name == genei) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

tiss_gt = tabix.read.table("../TFi-eQTL/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz",
                           paste0(strsplit(overlap_var,"_")[[1]][1],":",
                                  strsplit(overlap_var,"_")[[1]][2],"-",
                                  strsplit(overlap_var,"_")[[1]][2])) %>%
  select(-c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='indiv')

tiss_data = merge(tiss_tf_expr,tiss_gene_expr,
                  by=c('indiv','sample')) %>%
  merge(tiss_gt, by='indiv') %>%
  mutate(gt = factor(V1,
                     levels=c('0/0','0/1','1/1')))

tiss_corr = tiss_data %>%
  mutate(l10_v1y = log10(V1.y+0.001),
         l10_v1x = log10(V1.x+0.001)) %>%
  group_by(gt) %>%
  summarize(n_samp = n(),
            lm_b0 = lm(l10_v1y ~ l10_v1x, na.action = 'na.omit')$coefficients[1],
            lm_b1 = lm(l10_v1y ~ l10_v1x, na.action = 'na.omit')$coefficients[2])

tiss_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt)) +
  geom_abline(data=tiss_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','slateblue4','slateblue2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle(paste(sig_tiss_short, '\n', 
                descriptioni, genei, '\n', 
                overlap_var))

tiss_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=tiss_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','slateblue4','slateblue1')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle('Muscle eQTL')
ggsave(paste(paste0('plots/',tfi),descriptioni,overlap_var,"win_sig.pdf",sep='.'),
       height=2.5, width=3)




ns_tiss_short = 'ADPVSC'

nstiss_tf_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",ns_tiss_short,".tpm.gct"),
                          header=TRUE, sep='\t',
                          skip = 2) %>%
  filter(Description == tfi) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

nstiss_gene_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",ns_tiss_short,".tpm.gct"),
                            header=TRUE, sep='\t',
                            skip = 2) %>%
  filter(Name == genei) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

nstiss_gt = tabix.read.table("../TFi-eQTL/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz",
                           paste0(strsplit(overlap_var,"_")[[1]][1],":",
                                  strsplit(overlap_var,"_")[[1]][2],"-",
                                  strsplit(overlap_var,"_")[[1]][2])) %>%
  select(-c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='indiv')

nstiss_data = merge(nstiss_tf_expr,nstiss_gene_expr,
                  by=c('indiv','sample')) %>%
  merge(nstiss_gt, by='indiv') %>%
  mutate(gt = factor(V1,
                     levels=c('0/0','0/1','1/1')))

nstiss_corr = nstiss_data %>%
  mutate(l10_v1y = log10(V1.y+0.001),
         l10_v1x = log10(V1.x+0.001)) %>%
  group_by(gt) %>%
  summarize(n_samp = n(),
            lm_b0 = lm(l10_v1y ~ l10_v1x, na.action = 'na.omit')$coefficients[1],
            lm_b1 = lm(l10_v1y ~ l10_v1x, na.action = 'na.omit')$coefficients[2])

nstiss_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt)) +
  geom_abline(data=nstiss_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','chocolate4','chocolate2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle(paste(ns_tiss_short, '\n', 
                descriptioni, genei, '\n', 
                overlap_var))

nstiss_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=nstiss_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','chocolate4','chocolate2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle('Adipose eQTL')
ggsave(paste(paste0('plots/',tfi),descriptioni,overlap_var,"win_ns.pdf",sep='.'),
       height=2.5, width=3)
