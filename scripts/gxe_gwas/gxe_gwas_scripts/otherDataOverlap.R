#!/usr/bin/bash
##
##  otherDataOverlap.R
##
##  EDF 6/19/21
##

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggvenn)

setwd("~/projects/examples")

rep_data = read.table("~/data/Findley2021/cASE_GxE_rep.tab",
                      header=TRUE,sep='\t')

ct_hits = read.table("~/data/CTieQTLS/aaz8528_Kim-Hellmuth-Table-S1.csv",
                     skip=1, header=TRUE, sep=',') %>%
  unite(tiss_gene, tissue_id, gene_id, remove=FALSE)

sig_genes = read.table("replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                       header=TRUE, sep='\t') %>%
  separate(phenotype_id, c('ensg'), sep='[.]',
           remove=FALSE)
sig_data = read.table("replication/sig_assoc.fdr20.multi_or_cross_hits.vars_byset.txt",
                      header=TRUE, sep='\t') %>%
  separate(gene, c('ensg'), sep='[.]',
           remove=FALSE)

tiss_hits = read.table("replication/sig_assoc.fdr20.any_tiss_hits.all_hits.txt",
                       header=TRUE, sep='\t') %>%
  unite(tiss_gene, tiss, phenotype_id, remove=FALSE)

all_genes = read.table("variant_sets/caviar/caviar_var.uniqueeqtls.95set.txt",
                       header=TRUE,sep='\t')

all_genes_list = all_genes %>%
  select(ENSG00000000419.12) %>%
  unique() %>%
  separate(ENSG00000000419.12, c('ensg'), sep='[.]') %>%
  pull(ensg) %>%
  as.character()


## GxE stuff
table(unique(rep_data$ensg) %in% sig_genes$ensg)
table(unique(sig_genes$ensg) %in% rep_data$ensg)

table(all_genes_list %in% rep_data$ensg,
      all_genes_list %in% sig_genes$ensg)
table(all_genes_list %in% rep_data$ensg,
      all_genes_list %in% sig_genes$ensg) %>%
  fisher.test()


rep_sig_data = merge(rep_data,
      sig_genes,
      by='ensg')
nrow(rep_sig_data)
rep_sig_data %>%
  select(tf,Treatment,ensg) %>%
  unique() %>%
  group_by(tf,Treatment) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(-n)

genes_overlap_list = list(TF.eQTL = as.character(sig_genes$ensg),
                          GxE = as.character(rep_data$ensg))
ggvenn( genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('GxE_gene_venn.pdf',
       width=4, height=3)

rep_sig_data_table = rep_data %>% 
  select(ensg, SNP_Individual, Treatment) %>%
  unique() %>%
  merge(sig_data,
        by='ensg') %>%
  rename(GxE_SNP = SNP_Individual, GxE_Treatment=Treatment, top_TFeQTL_variant=variant) %>%
  select(gene, description, GxE_Treatment, GxE_SNP, tf, corr, top_TFeQTL_variant ) %>%
  arrange(gene, GxE_Treatment, tf)

write.table(rep_sig_data_table,
            "GxE_overlap.txt",
            col.names=TRUE,sep='\t',
            quote=FALSE, row.names=FALSE)




## CELL TYPE STUFF

table(unique(ct_hits$gene_id) %in% sig_genes$phenotype_id)
table(unique(sig_genes$phenotype_id) %in% ct_hits$gene_id)

table(unique(as.character(all_genes$ENSG00000000419.12)) %in% ct_hits$gene_id,
      unique(as.character(all_genes$ENSG00000000419.12)) %in% sig_genes$phenotype_id)
table(unique(as.character(all_genes$ENSG00000000419.12)) %in% ct_hits$gene_id,
      unique(as.character(all_genes$ENSG00000000419.12)) %in% sig_genes$phenotype_id) %>%
  fisher.test()

table(unique(ct_hits$gene_id) %in% tiss_hits$phenotype_id)
table(unique(tiss_hits$phenotype_id) %in% ct_hits$gene_id)
table(unique(ct_hits %>% 
               filter(tissue_id %in% as.character(tiss_hits$tiss)) %>% 
                        pull(tiss_gene)) %in%
             (tiss_hits %>%
                filter(tiss %in% as.character(ct_hits$tissue_id)) %>% 
                         pull(tiss_gene)))
table(unique(tiss_hits %>% 
               filter(tiss %in% as.character(ct_hits$tissue_id)) %>% 
               pull(tiss_gene)) %in%
        (ct_hits %>%
           filter(tissue_id %in% as.character(tiss_hits$tiss)) %>% 
           pull(tiss_gene)))
fisher.test(matrix(c(30000*9,4119,81,139), nrow=2))




tiss_hits %>%
  group_by(tiss) %>%
  summarize(n=n(),
            n_ct = sum(phenotype_id %in% 
                       (ct_hits %>% filter(tissue_id==as.character(unique(tiss))) %>% pull(gene_id) %>% as.character())),
            p_ct = n_ct/n,
            n_genes = length(unique(phenotype_id)),
            n_genes_ct = sum(unique(phenotype_id) %in% 
                               (ct_hits %>% filter(tissue_id==as.character(unique(tiss))) %>% pull(gene_id) %>% as.character())),
            p_genes_ct = n_genes_ct/n_genes)

tiss_hits %>%
  group_by(tiss, phenotype_id) %>%
  summarize(n_corrs = n()) %>%
  mutate(gene_ct = phenotype_id %in% ct_hits$gene_id && tiss %in% ct_hits$tiss) %>%
  ggplot(aes(gene_ct, n_corrs)) +
  geom_violin() +
  theme_classic()
tiss_hits %>%
  group_by(tiss, phenotype_id) %>%
  summarize(n_corrs = n()) %>%
  mutate(gene_ct = phenotype_id %in% ct_hits$gene_id && tiss %in% ct_hits$tiss) %>%
  ggplot(aes(gene_ct, n_corrs)) +
  geom_violin() +
  ylim(0,10) +
  theme_classic()
tiss_hits %>%
    group_by(tiss, phenotype_id) %>%
    summarize(n_corrs = n()) %>%
    mutate(gene_ct = phenotype_id %in% ct_hits$gene_id && tiss %in% ct_hits$tiss) %>%
    ggplot(aes(gene_ct, n_corrs)) +
    geom_violin() +
    ylim(10,NA) +
  theme_classic()

tiss_hits_ct = tiss_hits %>%
  group_by(tiss, phenotype_id) %>%
  summarize(n_corrs = n()) %>%
  mutate(gene_ct = phenotype_id %in% ct_hits$gene_id && tiss %in% ct_hits$tiss)
table(tiss_hits_ct$gene_ct)
wilcox.test(data = tiss_hits_ct, n_corrs ~ gene_ct)
tiss_hits_ct %>%
  group_by(gene_ct) %>%
  summarize(med=median(n_corrs),
            quart=quantile(n_corrs))

tiss_hits_ct_anytiss = tiss_hits %>%
  group_by(tiss, phenotype_id) %>%
  summarize(n_corrs = n()) %>%
  mutate(gene_ct = phenotype_id %in% ct_hits$gene_id)
table(tiss_hits_ct_anytiss$gene_ct)
wilcox.test(data = tiss_hits_ct_anytiss, n_corrs ~ gene_ct)
tiss_hits_ct_anytiss %>%
  group_by(gene_ct) %>%
  summarize(med=median(n_corrs),
            quart=quantile(n_corrs))


sum(table(ct_hits$tissue_id))
sum(table(ct_hits %>% 
            filter(tissue_id %in% c('Whole_Blood','Stomach','Colon_Transverse')) %>%
            pull(tissue_id)))




table(all_genes_list %in% rep_data$ensg,
      unique(as.character(all_genes$ENSG00000000419.12)) %in% ct_hits$gene_id)
table(all_genes_list %in% rep_data$ensg,
      unique(as.character(all_genes$ENSG00000000419.12)) %in% ct_hits$gene_id) %>%
  fisher.test()

