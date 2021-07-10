#!/usr/bin/Rscript
##
##  fig4_mocexamps.R
##
##  EDF 5/12/2021
##

library(UpSetR)
library(dplyr)
library(ggplot2)

setwd("~/projects/MANUSCRIPT/")

mpileups_comb_fi_test = read.table("data_tables/IRF1_kd/Fi_tst_combined_byvar_timecomb.allgenes.txt",
                                   header=TRUE, sep='\t')
mpileups_comb_fi_test_lowp = read.table("data_tables/IRF1_kd/Fi_tst_combined_byvar_timecomb.lowpgenes.txt",
                                        header=TRUE, sep='\t')

tested_genes = as.character(unique(mpileups_comb_fi_test$gene))
HEK_lowp_genes = as.character(unique(mpileups_comb_fi_test_lowp$gene))

IRF1_hits_cross = read.table("../overlap/replication/sig_assoc.fdr05.cross_hits.txt",
                             header=TRUE,sep='\t') %>%
  filter(phenotype_id %in% tested_genes,
         tf=='IRF1')

IRF1_hits_within = read.table("../overlap/replication/sig_assoc.fdr20.any_tiss_hits.txt",
                              header=TRUE, sep='\t' ) %>%
  filter(phenotype_id %in% tested_genes,
         tf=='IRF1')

IRF1_hits_multiorcross = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                                    header=TRUE, sep='\t') %>%
  filter(phenotype_id %in% tested_genes,
         tf=='IRF1')

pdf("plots/fig4_overlap.pdf",width=4,height=2)
upset(fromList(list(
  cross = IRF1_hits_cross %>% 
    pull(phenotype_id) %>% as.character(),
  tiss = IRF1_hits_within %>%
    pull(phenotype_id) %>% as.character(),
  HEK = HEK_lowp_genes)))
dev.off()

table(tested_genes %in% as.character(IRF1_hits_cross$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes))
fi_tst_cross = table(tested_genes %in% as.character(IRF1_hits_cross$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes)) %>% fisher.test

table(tested_genes %in% as.character(IRF1_hits_within$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes))
fi_tst_within = table(tested_genes %in% as.character(IRF1_hits_within$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes)) %>% fisher.test

table(tested_genes %in% as.character(IRF1_hits_multiorcross$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes))
table(tested_genes %in% as.character(IRF1_hits_multiorcross$phenotype_id),
      tested_genes %in% as.character(HEK_lowp_genes)) %>% fisher.test
fi_tst_moc = table(tested_genes %in% as.character(IRF1_hits_multiorcross$phenotype_id),
                   tested_genes %in% as.character(HEK_lowp_genes)) %>% fisher.test
data.frame(name=factor(c('Within', 'Cross', 'Dual Evidence'),
                       levels=c('Dual Evidence', 'Cross', 'Within')),
           OR = c(fi_tst_within$estimate, fi_tst_cross$estimate, fi_tst_moc$estimate),
                  logOR = c(log2(fi_tst_within$estimate), log2(fi_tst_cross$estimate), log2(fi_tst_moc$estimate)),
                  p = c(fi_tst_within$p.value, fi_tst_cross$p.value, fi_tst_moc$p.value)) %>%
  ggplot(aes(name,logOR)) +
  geom_col(width=.01) +
  geom_point() +
  geom_hline(yintercept=0) +
  theme_classic() +
  coord_flip()
ggsave("plots/fig4_overlap_enrich.pdf",
       height=1.5,
       width=2.5)

