#!/usr/bin/Rscript
##
##  calcOverlap_nouncorr.R
##
##
##  EDF 3/30/21
##

library(dplyr)
library(tidyr)


setwd("~/projects/overlap/")

args = commandArgs(trailingOnly=TRUE)

list = args[1]
#list="Thyroid.norm.ieqtl.sig_assoc.fdr20.txt"
#list="sig_assoc.fdr05.cross_hits.txt"
print(list)

# 
# tf_info = read.table("input_files/TFs.info.curated.txt",
#                      header=TRUE,sep='\t')

tf_files = list.files(path='perm_stats3',
                      pattern=paste0("*.perm_stats.",list),
                      full.names = TRUE)

all_overlap = do.call('rbind', lapply(tf_files, function(filei) {
  tf=strsplit(filei,split='[/.]')[[1]][2]
  datai = read.table(filei,
             header=TRUE,sep='\t') 
  
  if ( TRUE %in% datai$sig_gene ) {
    datai %>%
      pivot_longer(cols=c(eq2_chip,eq2_motif,eq2_both,
                          eq8_chip,eq8_motif,eq8_both)) %>%
      pivot_wider(id_cols=c(perm,name),
                  names_from=sig_gene,
                  values_from=c(value,num_genes)) %>%
      mutate(diff=value_TRUE-value_FALSE,
             true=value_TRUE,
             false=value_FALSE) %>%
      mutate(tf=tf)
  }
}))


overlap_perms = all_overlap %>%
  group_by(perm,name) %>%
  summarize(count_true = sum(!is.na(true)),
            count_false = sum(!is.na(false)),
            count_diff = sum(!is.na(diff)),
            true_stat = sum(true * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(true),0,num_genes_TRUE)),
            false_stat = sum(false * num_genes_FALSE, na.rm=TRUE) / sum(ifelse(is.na(false),0,num_genes_FALSE)),
            diff_stat = sum(diff * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(diff),0,num_genes_TRUE)),
            diff_stat_2 = true_stat - false_stat)


write.table(overlap_perms,
            paste0("sum_stats3/perms.",list),
            col.names=TRUE,row.names=FALSE,
            sep='\t',quote=FALSE)











# all_overlap = do.call('rbind', lapply(tf_files, function(filei) {
#   tf=strsplit(filei,split='[/.]')[[1]][2]
#   read.table(filei,
#              header=TRUE,sep='\t') %>%
#     mutate(tf=tf)
# }))
# 
# all_overlap %>% filter(sig_genes==FALSE)
# all_overlap %>% filter(perm==0)
# 
# sort(all_overlap %>% filter(perm==0, sig_genes) %>% pull(eq8_chip))
# sort(all_overlap %>% filter(perm==0, !(sig_genes)) %>% pull(eq8_chip))
# 
# sort(all_overlap %>% filter(perm==0) %>% pivot_wider(id_cols=c(perm,tf),names_from=sig_genes,values_from=eq8_chip) %>%
#        mutate(diff=`TRUE`-`FALSE`) %>% pull(diff))
# 
# sort(all_overlap %>% filter(perm==0, sig_genes) %>% pull(eq2_chip))
# sort(all_overlap %>% filter(perm==0, !(sig_genes)) %>% pull(eq2_chip))


# 
# 
# 
# 
# 
# test_files = c('sig_assoc.fdr20.multi_tiss_hits.txt',
#              'sig_assoc.fdr20.tiss_cross_hits.txt',
#              'sig_assoc.fdr20.multi_or_cross_hits.txt')
# tf_files2 = lapply(test_files, function(list) {
#   list.files(path='perm_stats',
#                       pattern=paste0("*.perm_stats.",list),
#                       full.names = TRUE)
# })
# 
# all_overlap2 = do.call('rbind', lapply(tf_files2, function(file_list) {
#   do.call('rbind', lapply(file_list, function(filei) {
#     tf=strsplit(filei,split='[/.]')[[1]][2]
#     list=strsplit(filei,split='[.]')[[1]][5]
#     read.table(filei,
#                header=TRUE,sep='\t') %>%
#       filter(sig_genes, perm==0) %>%
#       mutate(tf=tf,list=list)
#   }))
#   
# }))
# 
# 
# 
# 
