#!/usr/bin/Rscript
##
##  prot_stats.R
##
##
##  EDF 10/13/20
##

library(dplyr)
library(tidyr)
library(ggplot2)

tf_info = read.table("tf_measurements/TFs.info.curated.txt",
                     header=TRUE, sep='\t')

prot_corrs_stats = read.table("crosstiss_tf_corrs/correlations/cross_tiss_tf_protein_corrs.med0.MAF05.all_fdr20.stats.txt",
                              header=TRUE, sep='\t')
prot_vals = read.table("tf_measurements/protein_tissue_median",
                       header=TRUE, sep='\t')

prot_vals_tfs = prot_vals %>%
  filter(gene %in% (tf_info %>% separate(Name, c('Name'), '[.]') %>% pull(Name)))

prot_vals_tfs %>%
  pivot_longer(3:34) %>%
  filter(!is.na(value)) %>%
  group_by(c(as.character(description))) %>%
  summarise(n_tiss = n(),
            n_vals = length(unique(value))) %>%
  rename(tf=`c(as.character(description))`) %>%
  write.table("prot_stats.txt",
              col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)


all_prot_vals = read.table("prot_Jiang_2020/Table_S1_gene_in_protein_info.protein_normalized_abundance.csv",
                           header=TRUE,sep=',')
tf_genes = prot_vals_tfs %>% pull(gene)
all_prot_vals_summ = all_prot_vals %>%
  pivot_longer(cols=2:ncol(all_prot_vals)) %>%
  group_by(gene.id) %>%
  summarize(n=sum(!is.na(value)),
            min=min(value,na.rm=TRUE),
            max=max(value,na.rm=TRUE),
            var=var(value,na.rm=TRUE),
            sd=sqrt(var),
            med=median(value,na.rm=TRUE),
            mean=median(value,na.rm=TRUE),
            tf_gene = gene.id %in% as.character(tf_genes))
ggplot(all_prot_vals_summ,
       aes(tf_gene, log10(med))) +
  geom_boxplot() +
  theme_classic() +
  xlab('Protein is TF') +
  ylab('log10(median protein value across samples)')
all_prot_vals_summ %>%
  group_by(tf_gene) %>%
  summarize(median(med))
wilcox.test(med ~ tf_gene, all_prot_vals_summ)
wilcox.test(med ~ tf_gene, all_prot_vals_summ)$p.value

ggplot(all_prot_vals_summ,
       aes(tf_gene, log10(sd))) +
  geom_violin() +
  theme_classic()


expr_vals = read.table("tf_measurements/TFs.med_tpm.txt",
                       header=TRUE, sep='\t')
expr_vals %>%
  pivot_longer(cols=3:ncol(expr_vals)) %>%
  group_by(gene,description) %>%
  summarize(med=median(value),
            log2med = log2(med)) %>%
  ggplot(aes(log2med)) +
  geom_histogram()


