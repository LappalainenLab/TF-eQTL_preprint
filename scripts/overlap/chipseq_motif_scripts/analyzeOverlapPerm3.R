#!/usr/bin/Rscript
##
##  analyzeOvelapPerm.R
##
##  EDF 8/25/20
##

library(dplyr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(forcats)

setwd("~/projects/overlap/")

args = commandArgs(trailingOnly=TRUE)

gene_list = args[1]
print(gene_list)
# gene_list = "sig_assoc.fdr20.any_tiss_hits.txt"
# gene_list = "sig_assoc.fdr20.multi_tiss_hits.txt"
# gene_list = "sig_assoc.fdr20.tiss_cross_hits.txt"
# gene_list = "sig_assoc.fdr05.cross_hits.txt"
# gene_list = "sig_assoc.fdr20.protein_hits.txt"
#gene_list = "sig_assoc.fdr20.multi_or_cross_hits.txt"

tf_files = paste0("perm_stats3/",
                  list.files("perm_stats3/",
                             paste0("*.perm_stats.",gene_list)))

all_tf_data = do.call('rbind', lapply(tf_files, function(file_i) {
  tf_i=strsplit(basename(file_i),"[.]")[[1]][1]
  tf_data = read.table(file_i, header=TRUE, sep='\t') %>%
    mutate(tf=tf_i)
}) )





####################################### EQ 8 ###########################################


########### CHIP ###########################
eq8_chip = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq8_chip,num_genes)) %>%
  mutate(eq8_chip_diff=eq8_chip_TRUE-eq8_chip_FALSE,
         num_sig = num_genes_TRUE) 

eq8_chip_true = eq8_chip %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_chip_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq8_chip_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.chipseq.tfs.pdf"),
       height=10, width=14)
eq8_chip_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.chipseq.tfs.sortsig.pdf"),
       height=10, width=14)

eq8_chip_diff = eq8_chip %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_chip_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq8_chip_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.chipseq.tfs.pdf"),
       height=10, width=14)
eq8_chip_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.chipseq.tfs.sortsig.pdf"),
              height=10, width=14)


########### MOTIF ###########################
eq8_motif = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq8_motif,num_genes)) %>%
  mutate(eq8_motif_diff=eq8_motif_TRUE-eq8_motif_FALSE,
         num_sig = num_genes_TRUE) 

eq8_motif_true = eq8_motif %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_motif_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq8_motif_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.motif.tfs.pdf"),
       height=10, width=14)
eq8_motif_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.motif.tfs.sortsig.pdf"),
       height=10, width=14)

eq8_motif_diff = eq8_motif %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_motif_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq8_motif_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.motif.tfs.pdf"),
       height=10, width=14)
eq8_motif_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.motif.tfs.sortsig.pdf"),
       height=10, width=14)


########### BOTH ###########################
eq8_both = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq8_both,num_genes)) %>%
  mutate(eq8_both_diff=eq8_both_TRUE-eq8_both_FALSE,
         num_sig = num_genes_TRUE) 

eq8_both_true = eq8_both %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_both_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq8_both_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.both.tfs.pdf"),
       height=10, width=14)
eq8_both_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_true.both.tfs.sortsig.pdf"),
       height=10, width=14)

eq8_both_diff = eq8_both %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq8_both_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq8_both_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.both.tfs.pdf"),
       height=10, width=14)
eq8_both_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq8_diff.both.tfs.sortsig.pdf"),
       height=10, width=14)





######################## COMBINED ###########################

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq8_chip_true_all = sum(eq8_chip*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq8_chip_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("ChIPseq permutations\nAll TFs") +
  xlab("eq8 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq8_true.chipseq.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq8_chip_all = sum(eq8_chip*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq8_chip_all)) %>%
  mutate(eq8_chip_diff=`TRUE`-`FALSE`, 
         tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_chip_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`) %>%
  ggplot(aes(diff)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("ChIPseq permutations\nAll TFs") +
  xlab("eq8 diff stat")
ggsave(paste0("figs3/",gene_list,".eq8_diff.chipseq.all.pdf"),
       height=3, width=4)



all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq8_motif)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq8_motif_true_all = sum(eq8_motif*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq8_motif_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Motif permutations\nAll TFs") +
  xlab("eq8 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq8_true.motif.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq8_both)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq8_both_true_all = sum(eq8_both*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq8_both_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Both permutations\nAll TFs") +
  xlab("eq8 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq8_true.both.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq8_both)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq8_both_all = sum(eq8_both*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq8_both_all)) %>%
  mutate(eq8_both_diff=`TRUE`-`FALSE`, 
         tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_both_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`) %>%
  ggplot(aes(diff)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Both permutations\nAll TFs") +
  xlab("eq8 diff stat")
ggsave(paste0("figs3/",gene_list,".eq8_diff.both.all.pdf"),
       height=3, width=4)









####################################### EQ 2 ###########################################


########### CHIP ###########################
eq2_chip = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq2_chip,num_genes)) %>%
  mutate(eq2_chip_diff=eq2_chip_TRUE-eq2_chip_FALSE,
         num_sig = num_genes_TRUE) 

eq2_chip_true = eq2_chip %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_chip_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq2_chip_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.chipseq.tfs.pdf"),
       height=10, width=14)
eq2_chip_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.chipseq.tfs.sortsig.pdf"),
       height=10, width=14)

eq2_chip_diff = eq2_chip %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_chip_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq2_chip_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.chipseq.tfs.pdf"),
       height=10, width=14)
eq2_chip_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.chipseq.tfs.sortsig.pdf"),
       height=10, width=14)


########### MOTIF ###########################
eq2_motif = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq2_motif,num_genes)) %>%
  mutate(eq2_motif_diff=eq2_motif_TRUE-eq2_motif_FALSE,
         num_sig = num_genes_TRUE) 

eq2_motif_true = eq2_motif %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_motif_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq2_motif_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.motif.tfs.pdf"),
       height=10, width=14)
eq2_motif_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.motif.tfs.sortsig.pdf"),
       height=10, width=14)

eq2_motif_diff = eq2_motif %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_motif_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq2_motif_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.motif.tfs.pdf"),
       height=10, width=14)
eq2_motif_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.motif.tfs.sortsig.pdf"),
       height=10, width=14)


########### BOTH ###########################
eq2_both = all_tf_data %>%
  pivot_wider(id_cols = c(tf, perm), 
              names_from = sig_gene, 
              values_from = c(eq2_both,num_genes)) %>%
  mutate(eq2_both_diff=eq2_both_TRUE-eq2_both_FALSE,
         num_sig = num_genes_TRUE) 

eq2_both_true = eq2_both %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_both_TRUE) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`)


eq2_both_true %>% 
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.both.tfs.pdf"),
       height=10, width=14)
eq2_both_true %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(true)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_true.both.tfs.sortsig.pdf"),
       height=10, width=14)

eq2_both_diff = eq2_both %>%
  pivot_wider(id_cols = c(tf,num_sig), names_from = perm, values_from = eq2_both_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`)

eq2_both_diff %>% 
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~paste(tf,num_sig), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.both.tfs.pdf"),
       height=10, width=14)
eq2_both_diff %>% 
  as.data.frame() %>%
  arrange(-num_sig) %>%
  mutate(tf_order = factor(paste(as.character(tf), num_sig), 
                           levels=unique(paste(as.character(tf), num_sig))) ) %>%
  ggplot(aes(diff)) +
  geom_histogram(aes(fill=tf)) +
  geom_vline(aes(xintercept = real)) +
  facet_wrap(~tf_order, scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=5)) +
  theme(legend.position = 'none')
ggsave(paste0("figs3/",gene_list,".eq2_diff.both.tfs.sortsig.pdf"),
       height=10, width=14)





######################## COMBINED ###########################

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_chip_true_all = sum(eq2_chip*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_chip_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("ChIPseq permutations\nAll TFs") +
  xlab("eq2 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq2_true.chipseq.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_chip_all = sum(eq2_chip*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_chip_all)) %>%
  mutate(eq2_chip_diff=`TRUE`-`FALSE`, 
         tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = eq2_chip_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`) %>%
  ggplot(aes(diff)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("ChIPseq permutations\nAll TFs") +
  xlab("eq2 diff stat")
ggsave(paste0("figs3/",gene_list,".eq2_diff.chipseq.all.pdf"),
       height=3, width=4)



all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq2_motif)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_motif_true_all = sum(eq2_motif*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_motif_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Motif permutations\nAll TFs") +
  xlab("eq2 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq2_true.motif.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_motif_all = sum(eq2_motif*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_motif_all)) %>%
  mutate(eq2_motif_diff=`TRUE`-`FALSE`, 
         tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = eq2_motif_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`) %>%
  ggplot(aes(diff)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Motif permutations\nAll TFs") +
  xlab("eq2 diff stat")
ggsave(paste0("figs3/",gene_list,".eq2_diff.motif.all.pdf"),
       height=3, width=4)




all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq2_both)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_both_true_all = sum(eq2_both*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_both_true_all)) %>%
  mutate(tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = `TRUE`) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'true') %>%
  rename(real = `0`) %>%
  ggplot(aes(true)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Both permutations\nAll TFs") +
  xlab("eq2 corr genes stat")
ggsave(paste0("figs3/",gene_list,".eq2_true.both.all.pdf"),
       height=3, width=4)

all_tf_data %>%
  filter(num_genes > 0, num_genes < 32151,
         !is.na(eq2_both)) %>%
  group_by(perm, sig_gene) %>%
  summarize(eq2_both_all = sum(eq2_both*num_genes)/sum(num_genes)) %>%
  pivot_wider(id_cols = c(perm), 
              names_from = sig_gene, 
              values_from = c(eq2_both_all)) %>%
  mutate(eq2_both_diff=`TRUE`-`FALSE`, 
         tf='all') %>%
  pivot_wider(id_cols = tf, names_from = perm, values_from = eq2_both_diff) %>%
  pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
  rename(real = `0`) %>%
  ggplot(aes(diff)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle("Both permutations\nAll TFs") +
  xlab("eq2 diff stat")
ggsave(paste0("figs3/",gene_list,".eq2_diff.both.all.pdf"),
       height=3, width=4)







######################## Calc and plot perm p-vals EQ 8 ##########################



########### CHIP ############
eq8_chip_true_ps = eq8_chip_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_chip_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_chip_true_ps$p_perm_lt),-log10(eq8_chip_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_chip_true_ps$p_perm)
hist(eq8_chip_true_ps$p_perm_gt,breaks=20)
# 1-pi0est(eq8_chip_true_ps$p_perm)$pi0
# 1-pi0est(eq8_chip_true_ps$p_perm_lt)$pi0
# 1-pi0est(eq8_chip_true_ps$p_perm_gt)$pi0

write.table(eq8_chip_true_ps, "sum_stats3/multi_or_cross_hits.eq8_true.tf_chipseq_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)


eq8_chip_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_chip_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_chip_true_ps$num_sig,eq8_chip_true_ps$p_perm_gt,method='spearman')


eq8_chip_diff_ps = eq8_chip_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_chip_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_chip_diff_ps$p_perm_lt),-log10(eq8_chip_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_chip_diff_ps$p_perm)
hist(eq8_chip_diff_ps$p_perm_gt,breaks=20)
# 1-pi0est(eq8_chip_diff_ps$p_perm)$pi0
# 1-pi0est(eq8_chip_diff_ps$p_perm_lt)$pi0
# 1-pi0est(eq8_chip_diff_ps$p_perm_gt)$pi0

write.table(eq8_chip_diff_ps, "sum_stats3/multi_or_cross_hits.eq8_diff.tf_chipseq_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq8_chip_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_chip_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_chip_diff_ps$num_sig,eq8_chip_diff_ps$p_perm_gt,method='spearman')





########### MOTIF ############

eq8_motif_true_ps = eq8_motif_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_motif_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_motif_true_ps$p_perm_lt),-log10(eq8_motif_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_motif_true_ps$p_perm)
hist(eq8_motif_true_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq8_motif_true_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_motif_true_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_motif_true_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq8_motif_true_ps, "sum_stats3/multi_or_cross_hits.eq8_true.tf_motif_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq8_motif_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Motif permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_motif_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Motif permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_motif_true_ps$num_sig,eq8_motif_true_ps$p_perm_gt,method='spearman')


eq8_motif_diff_ps = eq8_motif_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_motif_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_motif_diff_ps$p_perm_lt),-log10(eq8_motif_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_motif_diff_ps$p_perm)
hist(eq8_motif_diff_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq8_motif_diff_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_motif_diff_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_motif_diff_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq8_motif_diff_ps, "sum_stats3/multi_or_cross_hits.eq8_diff.tf_motif_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq8_motif_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Motif permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_motif_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Motif permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_motif_diff_ps$num_sig,eq8_motif_diff_ps$p_perm_gt,method='spearman')




########### BOTH ############

eq8_both_true_ps = eq8_both_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_both_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_both_true_ps$p_perm_lt),-log10(eq8_both_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_both_true_ps$p_perm)
hist(eq8_both_true_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq8_both_true_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_both_true_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_both_true_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq8_both_true_ps, "sum_stats3/multi_or_cross_hits.eq8_true.tf_both_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq8_both_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Both permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_both_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Both permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_both_true_ps$num_sig,eq8_both_true_ps$p_perm_gt,method='spearman')


eq8_both_diff_ps = eq8_both_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq8_both_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq8_both_diff_ps$p_perm_lt),-log10(eq8_both_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq8_both_diff_ps$p_perm)
hist(eq8_both_diff_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq8_both_diff_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_both_diff_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq8_both_diff_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq8_both_diff_ps, "sum_stats3/multi_or_cross_hits.eq8_diff.tf_both_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq8_both_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Both permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq8_both_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Both permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq8_both_diff_ps$num_sig,eq8_both_diff_ps$p_perm_gt,method='spearman')





######################## Calc and plot perm p values: EQ 2 ##########################

########### CHIP ############
eq2_chip_true_ps = eq2_chip_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_chip_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_chip_true_ps$p_perm_lt),-log10(eq2_chip_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_chip_true_ps$p_perm)
hist(eq2_chip_true_ps$p_perm_gt,breaks=20)
# 1-pi0est(eq2_chip_true_ps$p_perm)$pi0
# 1-pi0est(eq2_chip_true_ps$p_perm_lt)$pi0
# 1-pi0est(eq2_chip_true_ps$p_perm_gt)$pi0

write.table(eq2_chip_true_ps, "sum_stats3/multi_or_cross_hits.eq2_true.tf_chipseq_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)


eq2_chip_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_chip_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_chip_true_ps$num_sig,eq2_chip_true_ps$p_perm_gt,method='spearman')


eq2_chip_diff_ps = eq2_chip_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_chip_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_chip_diff_ps$p_perm_lt),-log10(eq2_chip_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_chip_diff_ps$p_perm)
hist(eq2_chip_diff_ps$p_perm_gt,breaks=20)
# 1-pi0est(eq2_chip_diff_ps$p_perm)$pi0
# 1-pi0est(eq2_chip_diff_ps$p_perm_lt)$pi0
# 1-pi0est(eq2_chip_diff_ps$p_perm_gt)$pi0

write.table(eq2_chip_diff_ps, "sum_stats3/multi_or_cross_hits.eq2_diff.tf_chipseq_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq2_chip_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_chip_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-3,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("ChIPseq permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_chip_diff_ps$num_sig,eq2_chip_diff_ps$p_perm_gt,method='spearman')





########### MOTIF ############

eq2_motif_true_ps = eq2_motif_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_motif_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_motif_true_ps$p_perm_lt),-log10(eq2_motif_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_motif_true_ps$p_perm)
hist(eq2_motif_true_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq2_motif_true_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_motif_true_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_motif_true_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq2_motif_true_ps, "sum_stats3/multi_or_cross_hits.eq2_true.tf_motif_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq2_motif_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Motif permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_motif_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Motif permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_motif_true_ps$num_sig,eq2_motif_true_ps$p_perm_gt,method='spearman')


eq2_motif_diff_ps = eq2_motif_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_motif_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_motif_diff_ps$p_perm_lt),-log10(eq2_motif_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_motif_diff_ps$p_perm)
hist(eq2_motif_diff_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq2_motif_diff_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_motif_diff_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_motif_diff_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq2_motif_diff_ps, "sum_stats3/multi_or_cross_hits.eq2_diff.tf_motif_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq2_motif_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Motif permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_motif_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Motif permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_motif_diff_ps$num_sig,eq2_motif_diff_ps$p_perm_gt,method='spearman')




########### BOTH ############

eq2_both_true_ps = eq2_both_true %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(true <= real), ## how many perms are lt real
            num_gt = sum(true >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_both_true_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_both_true_ps$p_perm_lt),-log10(eq2_both_true_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_both_true_ps$p_perm)
hist(eq2_both_true_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq2_both_true_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_both_true_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_both_true_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq2_both_true_ps, "sum_stats3/multi_or_cross_hits.eq2_true.tf_both_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq2_both_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Both permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_both_true_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Both permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_both_true_ps$num_sig,eq2_both_true_ps$p_perm_gt,method='spearman')


eq2_both_diff_ps = eq2_both_diff %>%
  group_by(tf) %>%
  summarize(num_sig = first(num_sig),
            real_stat = first(real),
            num_lt = sum(diff <= real), ## how many perms are lt real
            num_gt = sum(diff >= real), ## how many perms are gt real
            p_perm = min( ( (min(num_lt,num_gt) + 1) / (10^4+1) ) * 2, 1),
            p_perm_lt = (num_lt + 1) / (10^4+1),
            p_perm_gt = (num_gt + 1) / (10^4+1))
qqplot(-log10(seq(2/(10^4+1),1,10^-4)),-log10(eq2_both_diff_ps$p_perm))
abline(a=0, b=1)
qqplot(-log10(eq2_both_diff_ps$p_perm_lt),-log10(eq2_both_diff_ps$p_perm_gt))
abline(a=0, b=1)
hist(eq2_both_diff_ps$p_perm)
hist(eq2_both_diff_ps$p_perm_gt,breaks=20)
tryCatch({1-pi0est(eq2_both_diff_ps$p_perm)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_both_diff_ps$p_perm_lt)$pi0},
         error=function(e) {print(e)} )
tryCatch({1-pi0est(eq2_both_diff_ps$p_perm_gt)$pi0},
         error=function(e) {print(e)} )

write.table(eq2_both_diff_ps, "sum_stats3/multi_or_cross_hits.eq2_diff.tf_both_stats.txt",
            col.names=TRUE,sep='\t',quote=FALSE,row.names=FALSE)

eq2_both_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [two-sided])') +
  ggtitle(paste("Both permutation p-values [two-sided] vs. number of corr. genes per TF\n",gene_list))
eq2_both_diff_ps %>%
  ggplot(aes(num_sig, -log10(p_perm_gt))) +
  geom_point(aes(col=tf)) +
  theme_classic() +
  geom_text(aes(label=ifelse(p_perm_gt<10^-2,tf,'')),
            hjust=0,vjust=0,size=2) +
  theme(legend.position='none') +
  xlab('Number correlated genes') +
  ylab('-log10(permutation p [one-sided])') +
  ggtitle(paste("Both permutation p-values [one-sided] vs. number of corr. genes per TF\n",gene_list))
cor.test(eq2_both_diff_ps$num_sig,eq2_both_diff_ps$p_perm_gt,method='spearman')










##################### Calc perm z scores.. not for now
# 
# 
# z_chip_all = all_tf_data %>%
#   filter(num_sig > 0) %>%
#   group_by(perm, sig_genes) %>%
#   summarize(eq8_chip_all = sum(eq8_chip*num_genes)/sum(num_genes)) %>%
#   pivot_wider(id_cols = c(perm), 
#               names_from = sig_genes, 
#               values_from = c(eq8_chip_all)) %>%
#   mutate(eq8_diff=`TRUE`-`FALSE`, 
#          tf='all') %>%
#   pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_diff) %>%
#   pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
#   rename(real = `0`) %>%
#   summarize(stat = 'eq8_chip_all',
#             list = gene_list,
#             mean=mean(diff),
#             sd = sd(diff),
#             real = first(real),
#             z = (first(real) - mean) / sd, 
#             num_lt = sum(real > diff),
#             num_gt = sum(diff > real),
#             p_perm = (min(num_lt,num_gt) + 1) / (10^4+1) * 2)
# z_motif_all = all_tf_data %>%
#   filter(num_sig > 0) %>%
#   filter(!is.na(eq8_motif)) %>%
#   group_by(perm, sig_genes) %>%
#   summarize(eq8_motif_all = sum(eq8_motif*num_genes)/sum(num_genes)) %>%
#   pivot_wider(id_cols = c(perm), 
#               names_from = sig_genes, 
#               values_from = c(eq8_motif_all)) %>%
#   mutate(eq8_diff=`TRUE`-`FALSE`, 
#          tf='all') %>%
#   pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_diff) %>%
#   pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
#   rename(real = `0`) %>%
#   summarize(stat = 'eq8_motif_all',
#             list = gene_list,
#             mean=mean(diff),
#             sd = sd(diff),
#             real = first(real),
#             z = (first(real) - mean) / sd,
#             num_lt = sum(real > diff),
#             num_gt = sum(diff > real),
#             p_perm = (min(num_lt,num_gt) + 1) / (10^4+1) * 2)
# z_both_all = all_tf_data %>%
#   filter(num_sig > 0) %>%
#   filter(!is.na(eq8_both)) %>%
#   group_by(perm, sig_genes) %>%
#   summarize(eq8_both_all = sum(eq8_both*num_genes)/sum(num_genes)) %>%
#   pivot_wider(id_cols = c(perm), 
#               names_from = sig_genes, 
#               values_from = c(eq8_both_all)) %>%
#   mutate(eq8_diff=`TRUE`-`FALSE`, 
#          tf='all') %>%
#   pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_diff) %>%
#   pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
#   rename(real = `0`) %>%
#   summarize(stat = 'eq8_both_all',
#             list = gene_list,
#             mean=mean(diff),
#             sd = sd(diff),
#             real = first(real),
#             z = (first(real) - mean) / sd,
#             num_lt = sum(real > diff),
#             num_gt = sum(diff > real),
#             p_perm = (min(num_lt,num_gt) + 1) / (10^4+1) * 2)
# z_all = rbind(z_chip_all, z_motif_all) %>%
#   rbind(z_both_all) 
# write.table(z_all, paste0("perm_stats3/z.",gene_list),
#               col.names=TRUE, row.names=FALSE,
#               sep='\t', quote=FALSE)
# 
# 
# 
# 
# tf_levels = eq8_chip_diff %>% arrange(real) %>% pull(tf) %>% as.character() %>% unique()
# eq8_chip_diff %>%
#   mutate(tf=factor(tf,levels=tf_levels)) %>%
#   group_by(tf,num_sig,real) %>%
#   filter(!is.na(real)) %>%
#   summarize(num_sig = first(num_sig),
#             num_lt = sum(diff < real), ## how many perms are lt real
#             num_gt = sum(diff > real), ## how many perms are gt real
#             p_perm = (1 - ( (max(num_lt,num_gt) ) / (10^4+1) ) ) * 2,
#             p_perm_lt = (num_lt + 1) / (10^4+1),
#             p_perm_gt = (num_gt + 1) / (10^4+1)) %>%
#   ggplot() +
#   geom_segment(aes(x=tf,xend=tf,y=0,yend=real)) +
#   geom_point(aes(x=tf,y=real,col=p_perm < 0.05,size=num_sig)) +
#   scale_color_manual(values=c('gray','darkorchid4')) +
#   theme_classic() +
#   coord_flip() +
#   ylab("Overlap enrichment statistic") +
#   ggtitle("ChIPseq overlap")
# ggsave(paste0("figs3/",gene_list,".chipseq.alltfs.pdf"),
#        height=8, width=4)
# 
# tf_levels = eq8_motif_diff %>% arrange(real) %>% pull(tf) %>% as.character() %>% unique()
# eq8_motif_diff %>%
#   mutate(tf=factor(tf,levels=tf_levels)) %>%
#   group_by(tf,num_sig,real) %>%
#   filter(!is.na(real)) %>%
#   summarize(num_sig = first(num_sig),
#             num_lt = sum(diff < real), ## how many perms are lt real
#             num_gt = sum(diff > real), ## how many perms are gt real
#             p_perm = (1 - ( (max(num_lt,num_gt) ) / (10^4+1) ) ) * 2,
#             p_perm_lt = (num_lt + 1) / (10^4+1),
#             p_perm_gt = (num_gt + 1) / (10^4+1)) %>%
#   ggplot() +
#   geom_segment(aes(x=tf,xend=tf,y=0,yend=real)) +
#   geom_point(aes(x=tf,y=real,col=p_perm < 0.05,size=num_sig)) +
#   scale_color_manual(values=c('gray','darkorchid4')) +
#   theme_classic() +
#   coord_flip() +
#   ylab("Overlap enrichment statistic") +
#   ggtitle("Motif overlap")
# ggsave(paste0("figs3/",gene_list,".motif.alltfs.pdf"),
#        height=8, width=4)
# 
# 
# tf_levels = eq8_both_diff %>% arrange(real) %>% pull(tf) %>% as.character() %>% unique()
# eq8_both_diff %>%
#   mutate(tf=factor(tf,levels=tf_levels)) %>%
#   group_by(tf,num_sig,real) %>%
#   filter(!is.na(real)) %>%
#   summarize(num_sig = first(num_sig),
#             num_lt = sum(diff < real), ## how many perms are lt real
#             num_gt = sum(diff > real), ## how many perms are gt real
#             p_perm = (1 - ( (max(num_lt,num_gt) ) / (10^4+1) ) ) * 2,
#             p_perm_lt = (num_lt + 1) / (10^4+1),
#             p_perm_gt = (num_gt + 1) / (10^4+1)) %>%
#   ggplot() +
#   geom_segment(aes(x=tf,xend=tf,y=0,yend=real)) +
#   geom_point(aes(x=tf,y=real,col=p_perm < 0.05,size=num_sig)) +
#   scale_color_manual(values=c('gray','darkorchid4')) +
#   theme_classic() +
#   coord_flip() +
#   ylab("Overlap enrichment statistic") +
#   ggtitle("Both overlap")
# ggsave(paste0("figs3/",gene_list,".both.alltfs.pdf"),
#        height=8, width=4)


# ### Investigate motif stuff
# 
# tf_motif_data = read.table("input_files/TFs.info.curated.txt",
#                            header=TRUE, sep='\t')
# tfs_a = tf_motif_data %>%
#   filter(Quality == "A") %>%
#   pull(TF)
# all_tf_data %>%
#   filter(tf %in% as.character(tfs_a)) %>%
#   group_by(perm, sig_genes) %>%
#   summarize(eq8_motif_all = sum(eq8_motif*num_genes, na.rm=TRUE) / sum(num_genes)) %>%
#   pivot_wider(id_cols = c(perm), 
#               names_from = sig_genes, 
#               values_from = c(eq8_motif_all)) %>%
#   mutate(eq8_diff=`TRUE`-`FALSE`, 
#          tf='all') %>%
#   pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_diff) %>%
#   pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
#   rename(real = `0`) %>%
#   ggplot(aes(diff)) +
#   geom_histogram() +
#   geom_vline(aes(xintercept = real)) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   ggtitle("Motif permutations -- \"A\" quality motif only")
# 
# all_tf_data %>%
#   filter(num_sig > 100) %>%
#   group_by(perm, sig_genes) %>%
#   summarize(eq8_motif_all = sum(eq8_motif*num_genes, na.rm=TRUE) / sum(num_genes)) %>%
#   pivot_wider(id_cols = c(perm), 
#               names_from = sig_genes, 
#               values_from = c(eq8_motif_all)) %>%
#   mutate(eq8_diff=`TRUE`-`FALSE`, 
#          tf='all') %>%
#   pivot_wider(id_cols = tf, names_from = perm, values_from = eq8_diff) %>%
#   pivot_longer(cols=as.character(1:10^4), names_to = c('perm'), values_to = 'diff') %>%
#   rename(real = `0`) %>%
#   ggplot(aes(diff)) +
#   geom_histogram() +
#   geom_vline(aes(xintercept = real)) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   ggtitle("Motif permutations\nAll TFs")
# 
# 
# eq8_motif_ps = eq8_motif_diff %>%
#   group_by(tf) %>%
#   summarize(num_lt = sum(real > diff),
#             num_gt = sum(diff > real),
#             p_perm = (min(num_lt,num_gt) + 1) / (10^4+1) * 2,
#             p_perm_lt = (num_lt + 1) / (10^4+1),
#             p_perm_gt = (num_gt + 1) / (10^4+1))
# qqplot(-log10(seq(2*10^-4,1,10^-4)),-log10(eq8_motif_ps$p_perm))
# abline(a=0, b=1)
# qqplot(-log10(seq(10^-4,1,10^-4)),-log10(eq8_motif_ps$p_perm_lt))
# abline(a=0, b=1)
# qqplot(-log10(seq(10^-4,1,10^-4)),-log10(eq8_motif_ps$p_perm_gt))
# abline(a=0, b=1)
# 
# print("DONE")

