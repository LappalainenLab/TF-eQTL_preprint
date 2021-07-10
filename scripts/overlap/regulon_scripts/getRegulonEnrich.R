#!/usr/bin/Rscript
##
##  getRegulonEnrich.R
##
##  EDF 4/12/2021
##

library(dplyr)
library(ggplot2)
library(tidyr)

setwd("~/projects/overlap/")


regulons = read.table("input_files/merged_interactions.txt",
                      header=TRUE,sep='\t')

regulons_hq = regulons %>%
  filter(database!='PAZAR') %>%
  group_by(TF,target) %>%
  summarize(cure = ifelse('curated_databases' %in% evidence, 1, 0),
            chip = ifelse('ChIP_Seq' %in% evidence, 1, 0),
            motif = ifelse('TFBS_scanning' %in% evidence, 1, 0),
            expr = ifelse('inferred' %in% evidence, 1, 0))

all_reg_genes = unique(regulons$target)

gene_info = read.table("input_files/gencode.v26.GRCh38.genes.gtf",
                       header=FALSE,sep='\t') %>%
  filter(V3=='gene') %>%
  separate(V9,sep=";",into=c('geneid','transcriptid','genetype','genename')) %>%
  separate(geneid,into=c(NA,'geneid'),sep=' ') %>%
  separate(genename,into=c(NA,NA,'genename'),sep=' ') %>%
  separate(genetype,into=c(NA,NA,'genetype'),sep=' ')

  


tf_genes_multiorcross = read.table("replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                                   header=TRUE,sep='\t') %>%
  merge(gene_info %>% select(geneid,genename,genetype),
        by.x='phenotype_id',
        by.y='geneid')


reg_fi_tsts = do.call('rbind', lapply(unique(as.character(tf_genes_multiorcross$tf)), function(tfi) {
  hitsi = filter(tf_genes_multiorcross,
                 tf == tfi)
  regsi = filter(regulons_hq,
                 TF == tfi)
  
  n_corr = nrow(hitsi)
  
  reg_any = regsi
  n_reg_any = nrow(reg_any)
  n_corr_reg_any = sum(as.character(hitsi$genename) %in% reg_any$target)
  fi_tst_any = fisher.test(table(factor(all_reg_genes %in% as.character(hitsi$genename),levels=c(TRUE,FALSE)),
                                 factor(all_reg_genes %in% as.character(reg_any$target), levels = c(TRUE,FALSE))))

  
  reg_cure = regsi %>% filter(cure==1)
  n_reg_cure = nrow(reg_cure)
  n_corr_reg_cure = sum(as.character(hitsi$genename) %in% reg_cure$target)
  fi_tst_cure = fisher.test(table(factor(all_reg_genes %in% as.character(hitsi$genename),levels=c(TRUE,FALSE)),
                                 factor(all_reg_genes %in% as.character(reg_cure$target), levels = c(TRUE,FALSE))))
  
  
  reg_hq = regsi %>% filter(cure==1 | chip==1)
  n_reg_hq = nrow(reg_hq)
  n_corr_reg_hq = sum(as.character(hitsi$genename) %in% reg_hq$target)
  fi_tst_hq = fisher.test(table(factor(all_reg_genes %in% as.character(hitsi$genename),levels=c(TRUE,FALSE)),
                                  factor(all_reg_genes %in% as.character(reg_hq$target), levels = c(TRUE,FALSE))))
  
  
  data.frame(tf=tfi,
             n_corr,
             reg=c('any','curated','cur_or_chip'),
             n_reg = c(n_reg_any,n_reg_cure,n_reg_hq),
             n_reg_corr = c(n_corr_reg_any,n_corr_reg_cure,n_corr_reg_hq),
             fi_OR = c(fi_tst_any$estimate,fi_tst_cure$estimate,fi_tst_hq$estimate),
             fi_logOR = log2(c(fi_tst_any$estimate,fi_tst_cure$estimate,fi_tst_hq$estimate)),
             fi_p = c(fi_tst_any$p.value,fi_tst_cure$p.value,fi_tst_hq$p.value))
} ) )


set.seed(9000)
reg_fi_tsts %>%
  filter(reg=='any') %>%
  ggplot(aes(1, fi_logOR)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=fi_p),
                 position=position_jitter(.2)) +
  theme_classic()

reg_fi_tsts %>% filter(reg=='any') %>% arrange(-fi_logOR) %>% head(10)



reg_fi_tsts %>%
  filter(reg=='cur_or_chip') %>%
  ggplot(aes(1, fi_logOR)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=fi_p),
             position=position_jitter(.2)) +
  theme_classic()

reg_fi_tsts %>% filter(reg=='cur_or_chip') %>% arrange(-fi_logOR) %>% head(10)




reg_fi_tsts %>%
  group_by(reg) %>%
  summarise(n_corr_tot = sum(n_corr),
            n_reg_tot = sum(n_reg),
            n_reg_corr_tot = sum(n_reg_corr),
            fi_OR_tot = fisher.test(
              matrix(c(length(all_reg_genes)*n() - n_corr_tot - n_reg_tot + n_reg_corr_tot,
                       n_corr_tot - n_reg_corr_tot,
                       n_reg_tot - n_reg_corr_tot,
                       n_reg_corr_tot), nrow = 2))$estimate,
            fi_p_tot = fisher.test(
              matrix(c(length(all_reg_genes)*n() - n_corr_tot - n_reg_tot + n_reg_corr_tot,
                       n_corr_tot - n_reg_corr_tot,
                       n_reg_tot - n_reg_corr_tot,
                       n_reg_corr_tot), nrow = 2))$p.value)



reg_fi_tsts_bycat = do.call('rbind', lapply(unique(as.character(tf_genes_multiorcross$tf)), function(tfi) {
  #print(tfi)
  hitsi = filter(tf_genes_multiorcross,
                 tf == tfi)
  regsi = filter(regulons_hq,
                 TF == tfi) %>%
    ungroup()
  
  n_corr = nrow(hitsi)
  
  if (nrow(regsi) > 0) 
  { 
    do.call('rbind', lapply(c('cure','chip','motif','expr'), function(defi) {
      #print(defi)
      
      reg_j = regsi %>%
        filter(select(regsi,matches(defi)) == 1)
      
      n_reg_j = nrow(reg_j)
      n_corr_reg_j = sum(as.character(hitsi$genename) %in% reg_j$target)
      fi_tst_j = fisher.test(table(factor(all_reg_genes %in% as.character(hitsi$genename), levels = c(TRUE,FALSE)),
                                   factor(all_reg_genes %in% as.character(reg_j$target), levels = c(TRUE,FALSE))))
      
      data.frame(tf=tfi,
                 n_corr,
                 reg=defi,
                 n_reg = n_reg_j,
                 n_reg_corr = n_corr_reg_j,
                 fi_OR = fi_tst_j$estimate,
                 fi_logOR = log2(fi_tst_j$estimate),
                 fi_p = fi_tst_j$p.value)
      
    }) )
  }
} ) )


reg_fi_tsts_overall = rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  group_by(reg) %>%
  summarize(n_tfs = n(),
            n_tot = n_tfs*length(all_reg_genes),
            n_corr = sum(n_corr),
            n_reg = sum(n_reg),
            n_reg_corr = sum(n_reg_corr),
            fi_OR = fisher.test(matrix(c(n_tot-n_corr-n_reg+n_reg_corr,
                                         n_corr-n_reg_corr,
                                         n_reg-n_reg_corr,
                                         n_reg_corr), nrow=2))$estimate,
            fi_log_OR = log2(fi_OR),
            fi_p = fisher.test(matrix(c(n_tot-n_corr-n_reg+n_reg_corr,
                                        n_corr-n_reg_corr,
                                        n_reg-n_reg_corr,
                                        n_reg_corr), nrow=2))$p.value) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr')))


reg_fi_tsts_bycat %>%
  group_by(reg) %>%
  summarise(n_corr_tot = sum(n_corr),
            n_reg_tot = sum(n_reg),
            n_reg_corr_tot = sum(n_reg_corr),
            fi_OR_tot = fisher.test(
              matrix(c(length(all_reg_genes)*n() - n_corr_tot - n_reg_tot + n_reg_corr_tot,
                       n_corr_tot - n_reg_corr_tot,
                       n_reg_tot - n_reg_corr_tot,
                       n_reg_corr_tot), nrow = 2))$estimate,
            fi_p_tot = fisher.test(
              matrix(c(length(all_reg_genes)*n() - n_corr_tot - n_reg_tot + n_reg_corr_tot,
                       n_corr_tot - n_reg_corr_tot,
                       n_reg_tot - n_reg_corr_tot,
                       n_reg_corr_tot), nrow = 2))$p.value)

reg_fi_tsts_bycat %>%
  ggplot(aes(reg,fi_logOR)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=fi_p),
             position=position_jitter(.2)) +
  theme_classic()

rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  mutate(fi_logOR_plot = ifelse(fi_OR == 0, -3, fi_logOR)) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  group_by(reg_lab) %>%
  mutate(fi_OR_avg = mean(fi_OR)) %>%
  ggplot(aes(reg_lab,fi_logOR_plot)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=fi_p,
                 shape=fi_p < 0.05),
             position=position_jitter(.2),
             size=.75) +
  scale_color_continuous(low='coral',high='gray50') +
  scale_shape_manual(values=c(21,19)) +
  geom_point(data=reg_fi_tsts_overall,
             aes(x=reg_lab,y=log2(fi_OR))) +
  theme_classic() +
  scale_y_continuous(labels=c('-Inf',-2,0,2,4,6),
                     breaks=c(-3,-2,0,2,4,6)) +
  xlab("Regulon Source") +
  ylab("log2(Odds Ratio)")
ggsave("figs3/regulons.orange.20210524.pdf",
       height=3,width=4)


rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  mutate(fi_logOR_plot = ifelse(fi_OR == 0, -3, fi_logOR)) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  ggplot(aes(reg_lab,fi_logOR_plot)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=fi_p),
             position=position_jitter(.2),
             size=.75) +
  scale_color_continuous(low='darkorchid4',high='gray50') +
  theme_classic() +
  scale_y_continuous(labels=c('-Inf',-2,0,2,4,6),
                     breaks=c(-3,-2,0,2,4,6)) +
  xlab("Regulon Source") +
  ylab("log2(Odds Ratio)")
ggsave("figs3/regulons.purple.20210524.pdf",
       height=3,width=4)



rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  mutate(fi_logOR_plot = ifelse(fi_OR == 0, -3, fi_logOR)) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  ggplot(aes(reg_lab,fi_logOR_plot)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=-log10(fi_p)),
             position=position_jitter(.2),
             size=.75) +
  geom_point(data=reg_fi_tsts_overall,
             aes(reg_lab,fi_log_OR),
             size=2) +
  scale_color_gradient2(low='gray',mid='darkorchid4',high='darkorchid4',
                        midpoint=3) +
  theme_classic() +
  scale_y_continuous(labels=c('-Inf',-2,0,2,4,6),
                     breaks=c(-3,-2,0,2,4,6)) +
  xlab("Regulon Source") +
  ylab("log2(Odds Ratio)")


rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  mutate(fi_logOR_plot = ifelse(fi_OR == 0, -3, fi_logOR)) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  filter(reg_lab!='curated') %>%
  ggplot(aes(reg_lab,fi_logOR_plot)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  geom_point(aes(col=-log10(fi_p)),
             position=position_jitter(.2),
             size=.75) +
  geom_point(data=reg_fi_tsts_overall %>% filter(reg_lab!='curated'),
             aes(reg_lab,fi_log_OR),
             size=2) +
  scale_color_gradient2(low='gray',mid='darkorchid4',high='darkorchid4',
                       midpoint=3) +
  theme_classic() +
  scale_y_continuous(labels=c('-Inf',-2,0,2,4,6),
                     breaks=c(-3,-2,0,2,4,6)) +
  xlab("Regulon Source") +
  ylab("log2(Odds Ratio)")

rbind(reg_fi_tsts %>% filter(reg=='any'),
      reg_fi_tsts_bycat) %>%
  mutate(fi_logOR_plot = ifelse(fi_OR == 0, -3, fi_logOR)) %>%
  mutate(reg_lab = factor(ifelse(reg=='cure','curated',as.character(reg)),
                          levels = c('any','curated','chip','motif','expr'))) %>%
  write.table("sum_stats3/regulon_enrich_bycat.txt",
              col.names=TRUE,row.names=FALSE,
              sep='\t',quote=FALSE)
write.table(reg_fi_tsts_overall,
            "sum_stats3/regulon_enrich_overall.txt",
            col.names=TRUE,row.names=FALSE,
            sep='\t',quote=FALSE)


