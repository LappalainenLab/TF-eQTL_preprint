#!/usr/bin/Rscript
##
##  combine_ADASTRAperms.R
##
##  ALT 01/31/2021
##

ADASTRA_tfs = read.table("ADASTRA/ADASTRA_tfs.txt", header=TRUE, sep="\t")

combined = as.data.frame(matrix(0, nrow=10001, ncol=7))
colnames(combined) = c("perm","sum_ASBstat_corr_05","sum_ASBstat_corr_10",
                       "num_corr_05", "num_corr_10", "avg_corr_05", "avg_corr_10")
combined$perm = c(0:10000)

for(tf in ADASTRA_tfs$TF){
  perms = read.table(paste0("ADASTRA_Susan/perm_stats_match/",tf,".perm_stats.txt"), header=TRUE, sep="\t")
  if(!is.na(perms$eq8_ASB_05[1]) & nrow(perms)==20002){
    for(i in c(1:10001)){
      perms_row_corr = perms[i*2,]
      num_corr = as.numeric(perms_row_corr$num_genes)
      ASBstat_corr_05 = perms_row_corr$eq8_ASB_05 * num_corr

      combined$sum_ASBstat_corr_05[i] = combined$sum_ASBstat_corr_05[i] + ASBstat_corr_05
      combined$num_corr_05[i] = combined$num_corr_05[i] + num_corr
    }
  }
  if(!is.na(perms$eq8_ASB_10[1]) & nrow(perms)==20002){
    for(i in c(1:10001)){
      perms_row_corr = perms[i*2,]
      num_corr = as.numeric(perms_row_corr$num_genes)
      ASBstat_corr_10 = perms_row_corr$eq8_ASB_10 * num_corr
      
      combined$sum_ASBstat_corr_10[i] = combined$sum_ASBstat_corr_10[i] + ASBstat_corr_10
      combined$num_corr_10[i] = combined$num_corr_10[i] + num_corr
    }
  }
}

combined$avg_corr_05 = combined$sum_ASBstat_corr_05 / combined$num_corr_05
combined$avg_corr_10 = combined$sum_ASBstat_corr_10 / combined$num_corr_10

write.table(combined,
            "ADASTRA_Susan/perm_stats_match/ADASTRA_Susan_combined_eq8.txt",
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## plot for combined across all ASB tfs
library(ggplot2)
min_05 = min(combined$avg_corr_05) - 0.005
max_05 = max(combined$avg_corr_05) + 0.005
scale_05 = (max_05 - min_05) / 100
plot = ggplot(combined, aes(avg_corr_05)) + geom_histogram(position="dodge",
                                                    breaks = seq(min_05, max_05, by=scale_05),
                                                    fill = "slategray3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("ADASTRA Susan (combined, , fdr05)") +
  geom_vline(xintercept = combined$avg_corr_05[1])
ggsave(paste0("ADASTRA_Susan/perms_combined_05_eq8.pdf"),
       width = 7.5, height = 5)
obs_avg_05 = combined$avg_corr_05[1]
tail_05 = nrow(combined[combined$avg_corr_05 >= obs_avg_05,])
p_value_05 = 2*tail_05 / 10001  # 0.12

