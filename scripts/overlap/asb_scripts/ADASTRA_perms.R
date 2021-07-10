#!/usr/bin/Rscript
##
##  ADASTRA_perms.R
##
##  ALT 01/19/2021
##

library(dplyr)
library(tidyr)

# args = commandArgs(trailingOnly = TRUE)
# tf = args[1]

setwd("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB")

# list of all accB variants
# accB = read.table(paste0("ADASTRA_Susan/by_tf/",tf,"_overview.txt"), header=TRUE, sep="\t", na.strings="NA")
accB = read.table(paste0("ADASTRA_Susan/combined_ASB/combined.txt"), header=TRUE, sep="\t", na.strings="NA")

# number of perms to run
perm_count = 10^4

# get total number of variants (unique var/gene pairs)
total_var = nrow(accB)

# calculate probability for any var to have ASB
expected_p_05 = nrow(accB[accB$has_ASB_fdr05==1,]) / total_var
expected_p_10 = nrow(accB[accB$has_ASB_fdr10==1,]) / total_var


print("Starting permutations...")
set.seed(90)
perm_info = do.call("rbind", lapply(1:perm_count, function(i) {
  if (i %% 1000 == 0){
    print(paste0("Now on perm ", i))
  }

  ## permute has_ASB; new column: has_ASB_perm
  perm = mutate(accB, has_ASB_perm_05 = sample(has_ASB_fdr05, n(), replace=FALSE),
                has_ASB_perm_10 = sample(has_ASB_fdr10, n(), replace=FALSE))

  ## calculate gene stats
  perm_stats = perm %>%
    group_by(gene) %>%
    summarize(n_var = n(),
              # tf = first(corr_tf),
              # tf_match = first(tf_match),
              corr_gene = first(corr_gene),
              p_ASB_05 = expected_p_05,
              obs_ASB_05 = sum(has_ASB_perm_05),
              exp_ASB_05 = expected_p_05 * n(),
              obs_ASB_10 = sum(has_ASB_perm_10),
              exp_ASB_10 = expected_p_10 * n(),
              diff_ASB_05 = (obs_ASB_05 - exp_ASB_05),
              diff_ASB_10 = (obs_ASB_10 - exp_ASB_10),
              stat2_ASB_05 = diff_ASB_05 / sqrt(exp_ASB_05),
              stat2_ASB_10 = diff_ASB_10 / sqrt(exp_ASB_10),
              perm = i ) %>%
    ungroup() %>%
    filter(n_var > 0) %>%
    as.data.frame()

  ## split genes into correlated/not and calculate means
  perm_stats %>%
    # group_by(perm, sig_genes = (tf_match==1)) %>%
    group_by(perm, sig_genes = (corr_gene==1)) %>%
    summarize(eq2_ASB_05 = mean(diff_ASB_05, na.rm=TRUE),
              eq2_ASB_10 = mean(diff_ASB_10, na.rm=TRUE),
              eq8_ASB_05 = mean(stat2_ASB_05, na.rm=TRUE),
              eq8_ASB_10 = mean(stat2_ASB_10, na.rm=TRUE),
              num_genes = n()) %>%
    group_by(perm) %>%
    mutate(num_sig = ifelse(n() == 2, min(num_genes), 0 ))
}) )


## calculate ASB stats for the original data, merge with perm data
all_info = accB %>%
  group_by(gene) %>%
  summarize(n_var = n(),
            # tf = first(corr_tf),
            # tf_match = first(tf_match),
            corr_gene = first(corr_gene),
            p_ASB_05 = expected_p_05,
            obs_ASB_05 = sum(has_ASB_fdr05),
            exp_ASB_05 = expected_p_05 * n(),
            obs_ASB_10 = sum(has_ASB_fdr10),
            exp_ASB_10 = expected_p_10 * n(),
            diff_ASB_05 = (obs_ASB_05 - exp_ASB_05),
            diff_ASB_10 = (obs_ASB_10 - exp_ASB_10),
            stat2_ASB_05 = diff_ASB_05 / sqrt(exp_ASB_05),
            stat2_ASB_10 = diff_ASB_10 / sqrt(exp_ASB_10),
            perm = 0 ) %>%
  ungroup() %>%
  filter(n_var > 0) %>%
  as.data.frame()  %>%
  # group_by(perm, sig_genes = (tf_match==1)) %>%
  group_by(perm, sig_genes = (corr_gene==1)) %>%
  summarize(eq2_ASB_05 = mean(diff_ASB_05, na.rm=TRUE),
            eq2_ASB_10 = mean(diff_ASB_10, na.rm=TRUE),
            eq8_ASB_05 = mean(stat2_ASB_05, na.rm=TRUE),
            eq8_ASB_10 = mean(stat2_ASB_10, na.rm=TRUE),
            num_genes = n()) %>%
  group_by(perm) %>%
  mutate(num_sig = ifelse(n() == 2, min(num_genes), 0 )) %>%
  rbind(perm_info)


## write all data to file
write.table(all_info,
            # paste0("ADASTRA_Susan/perm_stats_match/",tf,".perm_stats.txt"),
            "ADASTRA_Susan/combined_ASB/combined.perm_stats.txt",
            col.names = TRUE, row.names=FALSE,
            quote = FALSE, sep = '\t')



# ## plot perms
# library(ggplot2)
# library(gridExtra)
# 
# tf_list = read.table("TFs.info.curated.txt", header=TRUE, sep="\t")
# plots = lapply(tf_list$TF, function(tf){
#   try({
#     perms = read.table(paste0("ADASTRA_Susan/perm_stats_match/",tf,".perm_stats.txt"), header=TRUE,
#                        sep="\t", na.strings="NA")
# 
#     if(perms$num_sig[1] > 0){
#       perms_plot = perms[c(FALSE, TRUE),]
# 
#       min = min(perms_plot$eq2_ASB_05) - 0.05
#       max = max(perms_plot$eq2_ASB_05) + 0.05
#       scale = (max - min) / 100
#       return(ggplot(perms_plot, aes(eq2_ASB_05)) + geom_histogram(position="dodge",
#                                                             breaks = seq(min, max, by=scale),
#                                                             fill = "slategray3") +
#                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#                ggtitle(tf) + geom_vline(xintercept = perms_plot$eq2_ASB_05[1]))
#     }
#   })
# })
# 
# plots2 = plots[which(sapply(plots, is.list))]
# 
# pdf("ADASTRA_Susan/perms_plot_eq2_05.pdf", width = 40, height = 25)
# do.call(grid.arrange,plots2)
# dev.off()
# 
# 
# ## p-values
# ADASTRA_tfs = read.table("ADASTRA/ADASTRA_tfs.txt", header=TRUE, sep="\t")
# 
# p_values = data.frame(tf=character(), p_value=numeric(),
#                       num_sig=numeric(), obs_diff=numeric(), total_genes=numeric())
# 
# p_values = do.call("rbind", lapply(ADASTRA_tfs$TF, function(tf){
#   try({
#     perms = read.table(paste0("ADASTRA/perm_stats/",tf,".perm_stats.txt"), header=TRUE, sep="\t")
#     perms_diff = read.table(paste0("ADASTRA/perm_stats/",tf,".diff.txt"), header=TRUE, sep="\t")
#     obs_diff = perms_diff$diff[1]
#     num_sig = perms$num_sig[1]
#     total_genes = num_sig+perms$num_genes[1]
#     if(obs_diff >= median(perms_diff$diff)){
#       tail = nrow(perms_diff[perms_diff$diff >= obs_diff,])
#     }
#     else{
#       tail = nrow(perms_diff[perms_diff$diff <= obs_diff,])
#     }
#     p_value = 2*tail / 10001
#     return(data.frame(tf=tf, p_value=p_value, num_sig=num_sig,
#                       obs_diff=obs_diff, total_genes=total_genes))
#   })
# }))
# 
# p_values = p_values[!is.na(p_values$tf) & !is.na(p_values$obs_diff),]
# 
# ggplot(p_values, aes(x=num_sig, y=p_value)) + geom_point()
# ggplot(p_values, aes(x=obs_diff, y=p_value)) + geom_point()
# ggplot(p_values, aes(x=total_genes, y=p_value)) + geom_point()
# 
# 
# tfs = c("USF2", "MYC", "TFAP4", "NR4A1", "RFX1", "RELB", "ATF1")
# perms_plot = read.table(paste0("ADASTRA/perm_stats/",tf,".diff.txt"), header=TRUE,
#                         sep="\t", na.strings="NA")
# min = min(perms_plot$diff) - 0.05
# max = max(perms_plot$diff) + 0.05
# scale = (max - min) / 100
# ggplot(perms_plot, aes(diff)) + geom_histogram(position="dodge",
#                                                breaks = seq(min, max, by=scale),
#                                                fill = "slategray3") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   ggtitle(tf) + geom_vline(xintercept = perms_plot$diff[1])





# ## PLOTTING W CHIP-SEQ FILTER
# library(ggplot2)
# library(gridExtra)
# 
# tf_list = read.table("/gpfs/commons/groups/lappalainen_lab/eflynn/projects/overlap/sum_stats/multi_or_cross_hits.tf_chipseq_stats.txt",
#                      header=TRUE, sep="\t")
# tf_list = tf_list[!is.na(tf_list$real_stat),]
# tf_list = tf_list[tf_list$real_stat>0,]
# tf_list = tf_list[tf_list$p_perm > 0.05,]
# 
# plots = lapply(tf_list$tf, function(tf){
#   try({
#     perms = read.table(paste0("ADASTRA/perm_stats/",tf,".perm_stats.txt"), header=TRUE,
#                        sep="\t", na.strings="NA")
# 
#     if(perms$num_sig[1] > 0){
#       perms_plot = perms[c(FALSE, TRUE), c(1, 4)]
# 
#       min = min(perms_plot$eq8_ASB) - 0.05
#       max = max(perms_plot$eq8_ASB) + 0.05
#       scale = (max - min) / 100
#       return(ggplot(perms_plot, aes(eq8_ASB)) + geom_histogram(position="dodge",
#                                                             breaks = seq(min, max, by=scale),
#                                                             fill = "slategray3") +
#                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#                ggtitle(tf) + geom_vline(xintercept = perms_plot$eq8_ASB[1]))
#     }
#   })
# })
# 
# plots2 = plots[which(sapply(plots, is.list))]
# 
# pdf("ADASTRA/perms_plot_positive_chip_nonsig.pdf", width = 40, height = 25)
# do.call(grid.arrange,plots2)
# dev.off()
# 
# 
# combined = as.data.frame(matrix(0, nrow=10001, ncol=4))
# colnames(combined) = c("perm","sum_ASBstat_corr",
#                        "num_corr", "avg_corr")
# combined$perm = c(0:10000)
# 
# for(tf in tf_list$tf){
#   try({
#     perms = read.table(paste0("ADASTRA/perm_stats/",tf,".perm_stats.txt"), header=TRUE, sep="\t")
#     if(!is.na(perms$eq8_ASB[1]) & nrow(perms)==20002){
#       for(i in c(1:10001)){
#         perms_row_corr = perms[i*2,]
#         num_corr = as.numeric(perms_row_corr$num_genes)
#         ASBstat_corr = perms_row_corr$eq8_ASB * num_corr
#         
#         combined$sum_ASBstat_corr[i] = combined$sum_ASBstat_corr[i] + ASBstat_corr
#         combined$num_corr[i] = combined$num_corr[i] + num_corr
#       }
#     }
#   })
# }
# 
# combined$avg_corr = combined$sum_ASBstat_corr / combined$num_corr
# 
# # write.table(combined,
# #             "ADASTRA/perm_stats/ADASTRA_combined.txt",
# #             col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# 
# ## plot for combined across all ASB tfs
# min = min(combined$avg_corr) - 0.05
# max = max(combined$avg_corr) + 0.05
# scale = (max - min) / 100
# 
# obs_avg = combined$avg_corr[1]
# tail = nrow(combined[combined$avg_corr >= obs_avg,])
# p_value = 2*tail / 10001
# plot = ggplot(combined, aes(avg_corr)) + geom_histogram(position="dodge",
#                                                         breaks = seq(min, max, by=scale),
#                                                         fill = "slategray3") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   ggtitle("ADASTRA (chip: positive real stat, p_perm > 0.05)") +
#   labs(subtitle = paste0("p-value = ", p_value)) +
#   geom_vline(xintercept = combined$avg_corr[1])
# ggsave(paste0("ADASTRA/perms_combined_positive_chip_nonsig.pdf"),
#        width = 7.5, height = 5)
# 
# 
# tf_list = read.table("/gpfs/commons/groups/lappalainen_lab/eflynn/projects/overlap/sum_stats/multi_or_cross_hits.tf_chipseq_stats.txt",
#                      header=TRUE, sep="\t")
# tf_list = tf_list[!is.na(tf_list$real_stat),]
# 
# ASB_perms_info = do.call("rbind", lapply(tf_list$tf, function(tf){
#   try({
#     perms = read.table(paste0("ADASTRA/perm_stats/",tf,".perm_stats.txt"), header=TRUE,
#                        sep="\t", na.strings="NA")
#     perms = perms[perms$sig_genes,]
#     rownames(perms) = NULL
#     real_stat = perms$eq8_ASB[1]
#     right_tail = nrow(perms[perms$eq8_ASB >= real_stat,])
#     left_tail = nrow(perms[perms$eq8_ASB <= real_stat,])
#     if(right_tail < left_tail){
#       tail = right_tail
#     }
#     else{
#       tail = left_tail
#     }
#     p_value = 2*tail / 10001
#     if(p_value > 1){
#       p_value = 1
#     }
#     
#     return(data.frame(tf=tf, real_stat=real_stat, perm_p=p_value))
#   })
# }))
# 
# ASB_perms_info = ASB_perms_info[!is.na(ASB_perms_info$tf) & !is.na(ASB_perms_info$real_stat),]
# tf_list = tf_list[tf_list$tf %in% ASB_perms_info$tf,]
# 
# plot = cbind(tf_list, ASB_perms_info)
# plot = plot[,-c(9)]
# colnames(plot)[c(3,6,9,10)] = c("real_stat_chip", "p_perm_chip", "real_stat_ASB", "p_perm_ASB")
# plot$real_stat_ASB = as.numeric(plot$real_stat_ASB)
# plot$p_perm_ASB = as.numeric(plot$p_perm_ASB)
# 
# 
# library(ggplot2)
# ggplot(plot, aes(x=real_stat_chip, y=real_stat_ASB)) + geom_point()
# ggplot(plot, aes(x=p_perm_chip, y=p_perm_ASB)) + geom_point()
