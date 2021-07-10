#!/usr/bin/Rscript
##
##  ADASTRA_perms.R
##
##  ALT 03/08/2021
##

args = commandArgs(trailingOnly = TRUE)
num = args[1]
num = as.numeric(num)

setwd("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB")

ADASTRA_tfs = read.table("ADASTRA/ADASTRA_tfs.txt", header=TRUE, sep="\t")
accB_all = do.call("rbind", lapply(ADASTRA_tfs$TF, function(tf){
  return(read.table(paste0("ADASTRA_Susan/by_tf/",tf,"_overview.txt"), header=TRUE, sep="\t", na.strings="NA"))
}))

start = (500 * num) + 1
end = 500 * (num + 1)
if(end > length(unique(accB_all$gene))){
  end = length(unique(accB_all$gene))
}

accB = do.call("rbind", lapply(unique(accB_all$gene)[c(start:end)], function(gene){
  gene = as.character(gene)
  gene_rows = accB_all[accB_all$gene==gene,]
  add = do.call("rbind", lapply(unique(gene_rows$variant), function(var){
    var = as.character(var)
    has_ASB_fdr05 = 0
    has_ASB_fdr10 = 0
    corr_gene = 0
    rows = gene_rows[gene_rows$variant==var,]
    
    if(1 %in% rows$has_ASB_fdr05){  #ASB for any TF
      has_ASB_fdr05 = 1
    }
    if(1 %in% rows$has_ASB_fdr10){  #ASB for any TF
      has_ASB_fdr10 = 1
    }
    if(1 %in% rows$corr_gene){  #corr for any TF
      corr_gene = 1
    }
    corr_tf = rows$corr_tf[1]
    return(data.frame(variant = var, gene = gene, has_ASB_fdr05 = has_ASB_fdr05, 
                      has_ASB_fdr10 = has_ASB_fdr10, corr_gene = corr_gene))
  }))
  return(add)
}))

write.table(accB,
            paste0("ADASTRA_Susan/combined_ASB/combined_", num, ".txt"),
            col.names = TRUE, row.names=FALSE,
            quote = FALSE, sep = '\t')


combined = do.call("rbind", lapply(c(0:60), function(num){
  return(read.table(paste0("ADASTRA_Susan/combined_ASB/combined_", num, ".txt"), header=TRUE, sep="\t"))
}))
write.table(combined, "ADASTRA_Susan/combined_ASB/combined.txt", col.names = TRUE, row.names=FALSE,
            quote = FALSE, sep = '\t')

# plot perms (after running ADASTRA_perms.R)
library(ggplot2)

perms = read.table("ADASTRA_Susan/combined_ASB/combined.perm_stats.txt", header=TRUE,
                   sep="\t", na.strings="NA")
perms_plot = perms[c(FALSE, TRUE),] # only get data for corr genes
min = min(perms_plot$eq8_ASB_10) - 0.005
max = max(perms_plot$eq8_ASB_10) + 0.005
scale = (max - min) / 100
obs = perms_plot$eq8_ASB_10[1]
tail = nrow(perms_plot[perms_plot$eq8_ASB_10 >= obs,])
p_value = 2*tail / 10001
ggplot(perms_plot, aes(eq8_ASB_10)) + geom_histogram(position="dodge",
                                                  breaks = seq(min, max, by=scale),
                                                  fill = "slategray3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("ASB/corr any TF (eq8, 10% fdr) p-value=0.0002") + geom_vline(xintercept = perms_plot$eq8_ASB_10[1])
ggsave("ADASTRA_Susan/combined_ASB/perms_plot_10_eq8.pdf", width = 7.5, height = 5)


