#!/usr/bin/Rscript
##
##  calcOverlap3_comp.R
##
##
##  EDF 4/9/21
##

library(dplyr)
library(tidyr)


setwd("~/projects/overlap/")

args = commandArgs(trailingOnly=TRUE)

list = args[1]
#list="Thyroid.norm.ieqtl.sig_assoc.fdr20.txt"
#list="sig_assoc.fdr05.cross_hits.txt"
# list = "sig_assoc.fdr20.any_tiss_hits.txt"
print(list)

# 
tf_info = read.table("input_files/TFs.info.curated.txt",
                     header=TRUE,sep='\t')

tf_files = list.files(path='perm_stats3',
                      pattern=paste0("*.comp_stats.",list),
                      full.names = TRUE)

# tf_files = list.files(path='perm_stats3',
#                       pattern=paste0("*.comp_stats.txt"),
#                       full.names = TRUE)

all_overlap = do.call('rbind', lapply(tf_files, function(filei) {
  tf=strsplit(filei,split='[/.]')[[1]][2]
  datai = read.table(filei,
                     header=TRUE,sep='\t') 
  
  if ( TRUE %in% datai$sig_gene ) {
    datai %>%
      pivot_longer(cols=c(eq2_chip,eq2_motif,eq2_both,
                          eq8_chip,eq8_motif,eq8_both)) %>%
      pivot_wider(id_cols=c(tf,tf_comp,match,name),
                  names_from=sig_gene,
                  values_from=c(value,num_genes)) %>%
      mutate(diff=value_TRUE-value_FALSE,
             true=value_TRUE) %>%
      mutate(tf=tf)
  }
}))


overlap_comps = all_overlap %>%
  mutate(tf_match = ifelse(tf==as.character(tf_comp),'match',as.character(tf_comp))) %>%
  group_by(tf_match,name) %>%
  summarize(count_true = sum(!is.na(true)),
            count_diff = sum(!is.na(diff)),
            true_stat = sum(true * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(true),0,num_genes_TRUE)),
            diff_stat = sum(diff * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(diff),0,num_genes_TRUE)) )

# overlap_comps_long = overlap_comps %>%
#   pivot_longer(cols=c(true_stat,diff_stat),
#                names_to='var') %>%
#   pivot_wider(id_cols=c(name,var,count_true,count_diff),
#               names_from=tf_match) %>%
#   pivot_longer(cols=as.character(tf_info$TF),
#                names_to='tf_comp')
# 
# 
# overlap_comps_long %>%
#   ggplot() +
#   geom_hline(yintercept=0) +
#   geom_violin(aes(var,value)) +
#   geom_point(aes(var,match)) +
#   facet_wrap(~name) +
#   theme_classic()


write.table(overlap_comps,
            paste0("sum_stats3/comps.bytfcomp.",list),
            col.names=TRUE,row.names=FALSE,
            sep='\t',quote=FALSE)







overlap_comps_bytfexp = all_overlap %>%
  mutate(tf_match = ifelse(tf==as.character(tf_comp),'match',as.character(tf))) %>%
  group_by(tf,name) %>%
  summarize(count_true = sum(!is.na(true)),
            count_diff = sum(!is.na(diff)),
            true_stat = sum(true * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(true),0,num_genes_TRUE)),
            diff_stat = sum(diff * num_genes_TRUE, na.rm=TRUE) / sum(ifelse(is.na(diff),0,num_genes_TRUE)) )

write.table(overlap_comps,
            paste0("sum_stats3/comps.bytfexp.",list),
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


