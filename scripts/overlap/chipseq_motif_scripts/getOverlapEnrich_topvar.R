#!/usr/bin/Rscript
##
##  getOverlapEnrich_topvar.R 
##
##  EDF 6/23/21
##


library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/projects/overlap/")

args=commandArgs(trailingOnly = TRUE)

gene_list = args[1]
#gene_list = "replication/sig_assoc.fdr20.multi_or_cross_hits.vars_unique.txt"
print(gene_list)

## Table of TF name, etc.
tf_info = read.table("input_files/TFs.info.curated.txt",
                     header=TRUE,sep='\t')


sig_vars_temp = read.table(gene_list,
                          header=TRUE, sep='\t')
sig_vars_all = sig_vars_temp %>%
  unite(gene_var, 
        names(sig_vars_temp)[which(names(sig_vars_temp) %in% c('gene','phenotype_id'))],
        names(sig_vars_temp)[which(names(sig_vars_temp) %in% c('var','variant','var.min'))],
        remove=FALSE)


sig_var_overlap = do.call('rbind',apply(tf_info, 1, function(tf_row) {
  tf_i=as.character(tf_row[1])
  print(tf_i)
  
  sig_vars_tf = sig_vars_all %>%
    filter(tf==tf_i)
  
  tf_eqtls = read.table(paste0("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tf_i,".txt.gz"),
                         header=TRUE, sep='\t') %>%
    select(gene, gene_name, variant, tf, tf_chip, tf_motif) %>%
    unite(gene_var,gene,variant, remove=FALSE) %>%
    mutate(tf_both=tf_chip*tf_motif) %>%
    mutate(sig_var = gene_var %in% as.character(sig_vars_tf$gene_var)) %>%
    pivot_longer(cols=c(tf_chip,tf_motif,tf_both),
                 names_to='overlap')
  
  total_p_chip = tf_eqtls %>%
    filter(overlap=='tf_chip') %>%
    summarize( p = sum( value ) / 
                 n() ) %>%
    pull(p)
  total_p_motif = tf_eqtls %>%
    filter(overlap=='tf_motif') %>%
    summarize( p = sum( value ) /
                 n() ) %>%
    pull(p)
  total_p_both = tf_eqtls %>%
    filter(overlap=='tf_both') %>%
    summarize( p = sum( value ) /
                 n() ) %>%
    pull(p)

  tf_eqtls %>%
    group_by(sig_var, overlap) %>%
    summarize(.groups='drop',
              n_var = n(),
              tf = first(tf),
              p = ifelse(first(overlap)=='tf_chip', total_p_chip,
                         ifelse(first(overlap)=='tf_motif', total_p_motif,
                                ifelse(first(overlap)=='tf_both', total_p_both))),
              exp = p * n_var,
              obs = sum(value),
              p_obs = obs / n_var,
              diff = obs - exp,
              ratio = obs / exp)
}))

write.table(sig_var_overlap,
            paste0('top_var_overlap/overlap_stats.',strsplit(gene_list,"/")[[1]][2]),
            col.names=TRUE, row.names=FALSE,
            sep='\t',quote=FALSE)


sig_var_overlap %>%
  pivot_wider(id_cols = c(overlap, tf),
              names_from = sig_var,
              values_from = c(n_var, obs)) %>%
  group_by(overlap,tf) %>%
  mutate(fi_OR = ifelse(n_var_TRUE & obs_FALSE,
                        fisher.test(matrix(c(n_var_FALSE-obs_FALSE,
                               obs_FALSE,
                               n_var_TRUE-obs_TRUE,
                               obs_TRUE), nrow=2))$estimate,
                        NA),
         fi_logOR = ifelse(n_var_TRUE & obs_FALSE,
                           log2(fi_OR),
                           NA),
         fi_p = ifelse(n_var_TRUE & obs_FALSE,
                       fisher.test(matrix(c(n_var_FALSE-obs_FALSE,
                              obs_FALSE,
                              n_var_TRUE-obs_TRUE,
                              obs_TRUE), nrow=2))$p.value,
                       NA)) %>%
  write.table(paste0('top_var_overlap/fi_tsts_tf.',strsplit(gene_list,"/")[[1]][2]),
              col.names=TRUE, row.names=FALSE,
              sep='\t',quote=FALSE)

sig_var_overlap %>%
  pivot_wider(id_cols = c(overlap, tf),
              names_from = sig_var,
              values_from = c(n_var, obs)) %>%
  group_by(overlap) %>%
  summarize(n_var_FALSE_sum = sum(n_var_FALSE, na.rm=TRUE),
            n_var_TRUE_sum = sum(n_var_TRUE, na.rm=TRUE),
            obs_FALSE_sum = sum(obs_FALSE, na.rm=TRUE),
            obs_TRUE_sum = sum(obs_TRUE, na.rm=TRUE),
            fi_OR = fisher.test(matrix(c(n_var_FALSE_sum-obs_FALSE_sum,
                                         obs_FALSE_sum,
                                         n_var_TRUE_sum-obs_TRUE_sum,
                                         obs_TRUE_sum), nrow=2))$estimate,
            fi_logOR = log2(fi_OR),
            fi_p = fisher.test(matrix(c(n_var_FALSE_sum-obs_FALSE_sum,
                                        obs_FALSE_sum,
                                        n_var_TRUE_sum-obs_TRUE_sum,
                                        obs_TRUE_sum), nrow=2))$p.value) %>%
  write.table(paste0('top_var_overlap/fi_tsts_all.',strsplit(gene_list,"/")[[1]][2]),
              col.names=TRUE, row.names=FALSE,
              sep='\t',quote=FALSE)



