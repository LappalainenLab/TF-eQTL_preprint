#!/usr/bin/Rscript
##
##  examineTFCorrsMed0_byTF.R
##
##  EDF 5/11/20
##

library(dplyr)
library(ggplot2)
library(qvalue)
library(modelr)
library(tidyr)

setwd("~/projects/crosstiss_tf_corrs/")

args = commandArgs(trailingOnly = TRUE)
tf = args[1]
#tf = "ATF3"

print(paste("TF is",tf))


med0_tiss = read.table(paste0("correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tf,".txt.gz"),
                            header=TRUE, sep='\t')


hist(med0_tiss$num_tiss)
sum(!is.na(med0_tiss$sp_p))
ggplot(med0_tiss, aes(sp_p)) +
  geom_histogram() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position = 'none')


gene_meff = read.table("correlations/gene_tests.meff_gts_gao_99.txt",
                       header=TRUE, sep='\t')


med0_tiss_ps_meff = merge(med0_tiss, gene_meff, 
                          by='gene') %>%
  mutate( sp_p.bf_meff_gene = pmin(sp_p * meff, 1) )

# med0_tiss_ps_meff_all = med0_tiss_ps_meff %>%
#   group_by(tf) %>%
#   mutate( sp_p.bf_meff_gene.bh_tf = p.adjust(sp_p.bf_meff_gene, method='BH')) %>%
#   ungroup() %>%
#   mutate( sp_p.bf_meff_gene.bh_all = p.adjust(sp_p.bf_meff_gene, method='BH'))
# 
# table(med0_tiss_ps_meff_all$sp_p.bf_meff_gene.bh_all < 0.05)
# med0_tiss_ps_meff_all %>%
#   group_by(tf) %>%
#   summarize(sum(sp_p.bf_meff_gene.bh_tf < 0.05, na.rm=TRUE))
# med0_tiss_ps_meff_all %>%
#   group_by(tf) %>%
#   summarize(num_sig = sum(sp_p.bf_meff_gene.bh_tf < 0.05, na.rm=TRUE)) %>%
#   ggplot(aes(tf, num_sig)) +
#   geom_col(aes(fill=tf)) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ylim(c(0,8000))
# ggsave("plots/plot3.pdf")
#  
# write.table(med0_tiss_ps_meff_all, 
#              "correlations/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.adjp.txt.gz",
#              col.names=TRUE, row.names=FALSE,
#              sep='\t', quote=FALSE)

#################################################################################

med0_tiss_ps_meff_top = med0_tiss_ps_meff %>%
  group_by(gene) %>%
  filter(!is.na(sp_p)) %>%
  summarize(
    n_var = n(),
    meff = first(meff),
    sp_p.min = min(sp_p, na.rm=TRUE),
    sp_p.bf_meff_gene.min = min(sp_p.bf_meff_gene, na.rm=TRUE),
    sp_dir.min = ifelse(!is.na(sp_p.min), as.character(sp_dir[which.min(sp_p)]), NA),
    var.min = variant[which(sp_p == sp_p.min)][1],
    sp_rho.min = sp_rho[which(sp_p == sp_p.min)][1],
    tf_chip_sum = sum(tf_chip > 0, na.rm=TRUE),
    tf_motif_sum = sum(tf_motif > 0, na.rm=TRUE),
    tf_both_sum = sum(tf_chip > 0 & tf_motif > 0, na.rm=TRUE)
    ) %>%
  ungroup() %>%
  mutate( sp_p.bf_meff_gene.min.bh_tf = p.adjust(sp_p.bf_meff_gene.min, method='BH'),
          p_fdr_05_sig = sp_p.bf_meff_gene.min.bh_tf < 0.05 )

med0_tiss_ps_meff_top %>%
  ggplot(aes(sp_p.bf_meff_gene.min)) +
  geom_histogram() +
  theme_classic() +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle=45, hjust=1))


table(med0_tiss_ps_meff_top$sp_p.bf_meff_gene.min.bh_tf < 0.05)
table(med0_tiss_ps_meff_top$sp_p.bf_meff_gene.min.bh_tf < 0.05,
      med0_tiss_ps_meff_top$sp_dir.min)

write.table(med0_tiss_ps_meff_top,
            paste0("correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tf,".top.adjp.txt"),
            col.names=TRUE, row.names=FALSE,
            sep='\t', quote=FALSE)


tf_stats = med0_tiss_ps_meff_top %>%
  summarise(
    tf = tf,
    num = n(),
    num_sig = sum(p_fdr_05_sig),
    num_pos = sum(p_fdr_05_sig & sp_dir.min == 'pos'),
    num_neg = sum(p_fdr_05_sig & sp_dir.min == 'neg'),
    num_unc = sum(p_fdr_05_sig & sp_dir.min == 'unc'),
    pi1 = 1 - pi0est(sp_p.bf_meff_gene.min)$pi0.smooth[5],
    mean_num_var_sig = mean(ifelse(p_fdr_05_sig, n_var, NA), na.rm=TRUE),
    mean_num_var_ns = mean(ifelse(!(p_fdr_05_sig), n_var, NA), na.rm=TRUE),
    num_chip = sum(tf_chip_sum),
    num_motif = sum(tf_motif_sum),
    num_both = sum(tf_both_sum)
  )
write.table(tf_stats,
            paste0("correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tf,".stats.txt"),
            col.names=TRUE, row.names=FALSE,
            sep='\t', quote=FALSE)
