#!/usr/bin/Rscript
##
##  ADASTRA_overview.R
##
##  USAGE: ./ADASTRA_overview.R <tf>
##
##  ALT 01/13/2021
##

setwd("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB")

args = commandArgs(trailingOnly = TRUE)
tf = args[1]

# list of significant TF-gene pairs
sig_gene = read.table("../../../eflynn/projects/overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                      header=TRUE, sep="\t")

# list of all variants for each gene
all_variants = read.table("../../../eflynn/projects/TFi-eQTL/variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE, sep="\t")

# ADASTRA = read.delim(paste0("/gpfs/commons/groups/lappalainen_lab/data/ADASTRA/ADASTRA_TF_ford_dump/",tf,"_HUMAN.tsv"),
#                      header=TRUE, sep="\t", quote="")
ADASTRA = read.delim(paste0("/gpfs/commons/groups/lappalainen_lab/data/ADASTRA/ADASTRA_Susan/release_dump/TF/",tf,"_HUMAN.tsv"),
                     header=TRUE, sep="\t", quote="")

overview = do.call("rbind", lapply(unique(all_variants$gene), function(gene){
  gene = as.character(gene)

  add = data.frame(var=character(), gene=character(),
                   has_ASB_fdr05=numeric(), # 0 if variant doesn't have ASB; 1 if it does
                   has_ASB_fdr10=numeric(),
                   corr_gene=numeric(), # 0 if gene isn't correlated; 1 if correlated
                   corr_tf=character(), accB_tf=character(), tf_match = numeric(), # 0 if no match; 1 if match
                   fdrp_bh_ref=numeric(), fdrp_bh_alt=numeric())

  # get all variants for the gene
  variants = all_variants[all_variants$gene == gene,]

  # get all ADASTRA variants (accB) that match the gene's variants
  accB_var = do.call("rbind", lapply(variants$var, function(var){
    var = as.character(var)
    chr = strsplit(var, "_")[[1]][1]
    pos = strsplit(var, "_")[[1]][2]
    ref = strsplit(var, "_")[[1]][3]
    alt = strsplit(var, "_")[[1]][4]

    match = ADASTRA[ADASTRA$X.chr==chr & ADASTRA$pos==pos & ADASTRA$ref==ref & ADASTRA$alt==alt,]
    match = cbind(matrix(var, nrow=nrow(match), ncol=1), match)
    colnames(match)[1] = "variant"
    return(match)
  }))

  # if the gene has accB variant(s), add to dataframe
  if(nrow(accB_var) > 0){
    add = do.call("rbind", lapply(c(1:nrow(accB_var)), function(row){
      # automatically set tf_match and has_ASB to 0
      tf_match = 0
      has_ASB_fdr05 = 0
      has_ASB_fdr10 = 0

      accB_row = accB_var[row,]

      # check if variant has ASB
      if(accB_row$fdrp_bh_ref<0.05 | accB_row$fdrp_bh_alt<0.05){
        has_ASB_fdr05=1
      }
      if(accB_row$fdrp_bh_ref<0.1 | accB_row$fdrp_bh_alt<0.1){
        has_ASB_fdr10=1
      }

      # check if gene is correlated (in sig_gene) and get TF
      if(gene %in% sig_gene$phenotype_id){
        corr_gene=1

        # get all correlated TFs for the gene
        corr_tf = as.character(unique(sig_gene[sig_gene$phenotype_id==gene,]$tf))

        # check if ADASTRA and corr TF matches
        if(tf %in% corr_tf){
          tf_match=1
        }

        # combine corr_tf into a single string
        corr_tf = paste(corr_tf, sep="", collapse=" ")
      }
      else{
        corr_gene=0
        corr_tf="NA"
      }

      return(data.frame(variant=accB_row$variant, gene=gene, has_ASB_fdr05=has_ASB_fdr05,
                        has_ASB_fdr10=has_ASB_fdr10, corr_gene=corr_gene, corr_tf=corr_tf,
                        accB_tf=tf, tf_match=tf_match, fdrp_bh_ref=accB_row$fdrp_bh_ref,
                        fdrp_bh_alt=accB_row$fdrp_bh_alt))
    }))
  }
}))

write.table(overview, paste0("ADASTRA_Susan/by_tf/",tf,"_overview.txt"), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep="\t")




# ## plotting stats from ADASTRA overview
# library(ggplot2)
# ADASTRA = read.table("ADASTRA/ADASTRA_overview.txt", header=TRUE, sep="\t", na.strings="NA")
# ADASTRA_stats = data.frame(tf=character(), accB_var=numeric(), ASB_var=numeric(), accB_gene=numeric(),
#                            ASB_gene=numeric(), accB_corr_gene=numeric(), accB_corr_gene_tf_match=numeric(),
#                            ASB_corr_gene=numeric(), ASB_corr_gene_tf_match=numeric())
# ADASTRA_stats = do.call("rbind", lapply(unique(ADASTRA$accB_tf), function(tf){
#   tf = as.character(tf)
# 
#   ADASTRA_tf = ADASTRA[ADASTRA$accB_tf==tf,]
#   ADASTRA_tf_ASB = ADASTRA_tf[ADASTRA_tf$has_ASB==1,]
#   accB_var = nrow(ADASTRA_tf)
#   ASB_var = nrow(ADASTRA_tf_ASB)
#   accB_gene = length(unique(ADASTRA_tf$gene))
#   ASB_gene = length(unique(ADASTRA_tf_ASB$gene))
# 
#   ADASTRA_tf_corr = ADASTRA_tf[ADASTRA_tf$corr_gene==1,]
#   ADASTRA_tf_ASB_corr = ADASTRA_tf_ASB[ADASTRA_tf_ASB$corr_gene==1,]
#   accB_corr_gene = length(unique(ADASTRA_tf_corr$gene))
#   ASB_corr_gene = length(unique(ADASTRA_tf_ASB_corr$gene))
# 
#   ADASTRA_tf_corr_match = ADASTRA_tf_corr[ADASTRA_tf_corr$tf_match==1,]
#   ADASTRA_tf_ASB_corr_match = ADASTRA_tf_ASB_corr[ADASTRA_tf_ASB_corr$tf_match==1,]
#   accB_corr_gene_tf_match = length(unique(ADASTRA_tf_corr_match$gene))
#   ASB_corr_gene_tf_match = length(unique(ADASTRA_tf_ASB_corr_match$gene))
# 
#   return(data.frame(tf=tf, accB_var=accB_var, ASB_var=ASB_var, accB_gene=accB_gene,
#                     ASB_gene=ASB_gene, accB_corr_gene=accB_corr_gene,
#                     accB_corr_gene_tf_match=accB_corr_gene_tf_match, ASB_corr_gene=ASB_corr_gene,
#                     ASB_corr_gene_tf_match=ASB_corr_gene_tf_match))
# }))
# 
# write.table(ADASTRA_stats, "ADASTRA/ADASTRA_stats.txt", col.names=TRUE,
#             row.names=FALSE, quote=FALSE, sep="\t")
# 
# ggplot(ADASTRA_stats, aes(x=accB_gene, y=ASB_gene)) + geom_point()s



# p_values = read.table("ADASTRA/perm_stats/ADASTRA_p_values.txt", header=TRUE, sep="\t")
# min = -0.01
# max = 1.01
# scale = 0.05
# ggplot(p_values, aes(p_value)) + geom_histogram(position="dodge",
#                                                            breaks = seq(min, max, by=scale),
#                                                            fill = "slategray3") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   ggtitle("ADASTRA perm p-values")
# 
# ADASTRA_stats = ADASTRA_stats[ADASTRA_stats$tf %in% p_values$tf,]
# ADASTRA_stats = cbind(ADASTRA_stats, p_values$p_value, p_values$obs_diff)
# colnames(ADASTRA_stats)[c(11,12)] = c("p_value","obs_diff")
# colnames(ADASTRA_stats)[9] = "ASB_TFeQTL"
# 
# ADASTRA_stats %>%
#   mutate(tf = fct_reorder(tf, obs_diff)) %>%
#   filter(!is.na(obs_diff)) %>%
#   ggplot(aes(tf,obs_diff)) +
#   geom_col(width=.1) +
#   geom_point(aes(col=p_value < 0.05,size=ASB_TFeQTL)) +
#   theme_classic() +
#   scale_color_manual(values=c('gray','darkorchid4')) +
#   coord_flip()
# 
# ggsave("ADASTRA/ADASTRA_overview.pdf", height=12, width=6, units="in")
# 
# 
# ggplot(ADASTRA_stats, aes(x=exp_ASB,y=ASB_TFeQTL)) + ylab("obs_ASB") +
#   geom_point(aes(col=p_value))
# 
# ggplot(ADASTRA_stats, aes(x=exp_ASB,y=ASB_TFeQTL)) + ylab("obs_ASB") +
#   geom_point() + geom_abline(intercept = 0, slope = 1)
# 
# ggplot(ADASTRA_stats, aes(x=log10(exp_ASB),y=log10(ASB_TFeQTL))) + xlab("log10(exp_ASB)") +
#   ylab("log10(obs_ASB)") + xlim(-2.5, NA) + ylim(-2.5, NA) +
#   geom_point() + geom_abline(intercept = 0, slope = 1)
