#!/usr/bin/Rscript
##
##
##  permGenesTFs_bytf3.R
##
##  EDF 8/21/21
##

library(dplyr)
library(tidyr)

setwd("~/projects/overlap/")

args=commandArgs(trailingOnly = TRUE)

## List of gene/TF pairs that have significant correlations
gene_list = args[1]
# gene_list = "replication/sig_assoc.fdr20.any_tiss_hits.txt"

## TF of interest
tf_i = args[2]
# tf_i = "BATF"

## Number of permutations to run
perm_n = 10^4
# perm_n = 100

## Read gene list and filter only for genes correlated with tf of interest
sig_genes = read.table(gene_list,
                       header=TRUE, sep='\t') %>%
  filter(tf == tf_i)

## File with gene/variant/overlap info
all_eqtls = read.table(paste0("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tf_i,".txt.gz"),
                       header=TRUE, sep='\t') %>%
  select(gene, gene_name, variant, tf, tf_chip, tf_motif) %>%
  unite(tf_chip_motif, tf_chip, tf_motif, remove=FALSE) %>%
  mutate(sig_gene = gene %in% as.character(sig_genes$phenotype_id))

## Table of TF name, etc.
tf_info = read.table("input_files/TFs.info.curated.txt",
                     header=TRUE,sep='\t')


## 1. Calculate probability of TF overlaps:
total_p_chip = all_eqtls %>%
  summarize( p = sum( tf_chip > 0) / 
               n() ) %>%
  pull(p)

total_p_motif = all_eqtls %>%
  summarize( p = sum( tf_motif > 0) / 
               n() ) %>%
  pull(p)


total_p_both = all_eqtls %>%
  summarize( p = sum( tf_chip > 0 & tf_motif > 0 ) / 
               n() ) %>%
  pull(p)


## 2. Permute TF overlap & calc stats
print("Starting permutations...")
set.seed(90)
perm_info = do.call('rbind', lapply(1:perm_n, function(i) {
  if (i %% 1000 == 0) { print(paste("Now on perm",i)) }
  
  ## 2a. Permute TF overlap
  perm_i = all_eqtls %>% 
    mutate(tf_chip_motif_perm = sample(tf_chip_motif, replace=FALSE)) %>%
    mutate(tf_chip_perm = ifelse(startsWith(tf_chip_motif_perm,"1"), 1, 0),
           tf_motif_perm = ifelse(endsWith(tf_chip_motif_perm,"1"), 1, 0))
  # separate(tf_chip_motif_perm, into = c("tf_chip_perm", "tf_motif_perm"), sep='_' ) %>%
  # mutate(tf_chip_perm = as.numeric(tf_chip_perm), 
  #        tf_motif_perm = as.numeric(tf_motif_perm))
  
  ## 2b. Calculate gene stats
  perm_i_annot = perm_i %>%
    group_by(gene) %>%
    summarize(.groups = 'drop', 
              sig_gene = first(sig_gene),
              n_var = n(),
              #tf = tfi,
              #p_chip = total_p_chip,
              obs_chip = sum(tf_chip_perm),
              exp_chip = total_p_chip * n(),
              diff_chip = (obs_chip - exp_chip), 
              stat2_chip = diff_chip / sqrt(exp_chip),
              #p_motif = total_p_motif,
              obs_motif = sum(tf_motif_perm),
              exp_motif = total_p_motif * n(),
              diff_motif = (obs_motif - exp_motif),
              stat2_motif = diff_motif / sqrt(exp_motif),
              #p_both = total_p_both,
              obs_both = sum(tf_chip_perm & tf_motif_perm),
              exp_both = total_p_both * n(),
              diff_both = (obs_both - exp_both),
              stat2_both = diff_both / sqrt(exp_both),
              perm = i ) %>%
    #ungroup() %>%
    filter(n_var > 0)
  
  ## 2c. Split genes into correlated/not and calculate means
  perm_i_annot %>%
    group_by(perm, sig_gene) %>%
    summarize(.groups='drop',
              eq2_chip = mean(diff_chip, na.rm=TRUE),
              eq8_chip = mean(stat2_chip, na.rm=TRUE),
              eq2_motif = mean(diff_motif, na.rm=TRUE),
              eq8_motif = mean(stat2_motif, na.rm=TRUE),
              eq2_both = mean(diff_both, na.rm=TRUE),
              eq8_both = mean(stat2_both, na.rm=TRUE),
              num_genes = n()) #%>%
  #group_by(perm) %>%
  #mutate(num_sig = ifelse(n() == 2, min(num_genes), 0 ))
}) )


## 3. Calculate TF overlap stats for the original data, merge with permut. data
all_info = all_eqtls %>%
  group_by(gene) %>%
  summarize(.groups='drop',
            sig_gene = first(sig_gene),
            n_var = n(),
            tf = first(tf),
            p_chip = total_p_chip,
            obs_chip = sum(tf_chip),
            exp_chip = total_p_chip * n(),
            diff_chip = (obs_chip - exp_chip), 
            stat2_chip = diff_chip / sqrt(exp_chip),
            p_motif = total_p_motif,
            obs_motif = sum(tf_motif),
            exp_motif = total_p_motif * n(),
            diff_motif = (obs_motif - exp_motif),
            stat2_motif = diff_motif / sqrt(exp_motif),
            p_both = total_p_both,
            obs_both = sum(tf_chip & tf_motif),
            exp_both = total_p_both * n(),
            diff_both = (obs_both - exp_both),
            stat2_both = diff_both / sqrt(exp_both),
            perm = 0 ) %>%
  filter(n_var > 0) %>%
  group_by(perm, sig_gene) %>%
  summarize(.groups='drop',
            eq2_chip = mean(diff_chip, na.rm=TRUE),
            eq8_chip = mean(stat2_chip, na.rm=TRUE),
            eq2_motif = mean(diff_motif, na.rm=TRUE),
            eq8_motif = mean(stat2_motif, na.rm=TRUE),
            eq2_both = mean(diff_both, na.rm=TRUE),
            eq8_both = mean(stat2_both, na.rm=TRUE),
            num_genes = n()) %>%
  #group_by(perm) %>%
  #mutate(num_sig = ifelse(n() == 2, min(num_genes), 0 )) %>%
  rbind(perm_info)


## 4. Write all data to file
write.table(all_info, 
            paste0("perm_stats3/",tf_i,".perm_stats.",basename(gene_list)),
            col.names = TRUE, row.names=FALSE,
            quote = FALSE, sep = '\t')
