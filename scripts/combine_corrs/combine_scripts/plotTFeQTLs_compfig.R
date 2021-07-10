#!/usr/bin/Rscript
##
##  plotTFeQTLs_compfig.R
##
##  EDF 5/21/21
##

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(heatmap3)

setwd("~/projects/TFeQTLs_figs/")

tf_info = read.table("../TFi-eQTL/input_files/TFs.info.curated.txt",
                     header=TRUE, sep='\t')
tiss_info = read.table("../TFi-eQTL/input_files/tissue_info.final.txt",
                       header=TRUE, sep='\t')
genes = read.table("../TFi-eQTL/variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                   header=TRUE, sep='\t') %>%
  group_by(gene) %>%
  summarize(num=n())


sig_hits = read.table("../TFi-eQTL/eqtl_sig/all.fdr20.sig_eqtls.txt",
                      header=TRUE, sep='\t') %>%
  merge(tf_info %>% select(1:2) %>% rename(tf_id = Name),
        by.x=c('tf'), by.y=c('TF')) %>%
  merge(tiss_info, by.x="tiss", by.y="TISSUE_NAME") %>%
  filter(tf_id != as.character(gene)) %>%
  mutate(tiss = factor(tiss, levels=unique(as.character(tiss_info$TISSUE_NAME))))


all_genes = genes %>% pull(gene) %>% as.character()
all_tfs = tf_info %>% pull(TF) %>% as.character()
all_gene_tf = sapply(all_genes, function(gene) {
  sapply(all_tfs, function(tf) {
    paste(gene,tf,sep='_')
  })
})

fi_tsts_fdr20 = do.call('rbind', 
                        lapply(unique(as.character(sig_hits$TISSUE_ABBRV)), function(tiss1) {
                          do.call('rbind',
                                  lapply(unique(as.character(sig_hits$TISSUE_ABBRV)), function(tiss2) {
                                    print(c(tiss1,tiss2))
                                    this_table=table(all_gene_tf %in% (sig_hits %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                                     all_gene_tf %in% (sig_hits %>% filter(TISSUE_ABBRV==tiss2) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)))
                                    this_fi = fisher.test(this_table)
                                    print(this_fi)
                                    data.frame(tiss1, tiss2,
                                               fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                                  }) )
                        } ) )


fi_tsts_fdr20_mat = fi_tsts_fdr20 %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_fdr20_mat_mod = apply(fi_tsts_fdr20_mat, c(1,2), function(x) {
  max(min(x, 14), 0)
})
heatmap3(fi_tsts_fdr20_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2')



sig_hits_cross = read.table("../overlap/replication/sig_assoc.fdr05.cross_hits.txt",
                            header=TRUE, sep='\t')
fi_tsts_cross = do.call('rbind', 
                        lapply(unique(as.character(sig_hits$TISSUE_ABBRV)), function(tiss1) {
                          print(c(tiss1,'cross'))
                          this_table=table(all_gene_tf %in% (sig_hits %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                           all_gene_tf %in% (sig_hits_cross %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)))
                          this_fi = fisher.test(this_table)
                          print(this_fi)
                          data.frame(tiss1=c(tiss1), tiss2=c('cross'),
                                     fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                        } ) )
fi_tsts_cross_mat = fi_tsts_cross %>%
  filter(tiss2=='cross') %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_cross_mat_mod = apply(fi_tsts_cross_mat, c(1,2), function(x) {
  max(min(x, 13), -4)
})

this_table=table(all_gene_tf %in% (sig_hits %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                 all_gene_tf %in% (sig_hits_cross %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)))
this_fi = fisher.test(this_table)
print(this_fi)
data.frame(tiss1=c('all_tiss'), tiss2=c('cross'),
           fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
.Machine$double.xmin



sig_hits_prot = read.table("../overlap/replication/sig_assoc.fdr20.protein_hits_filtered.txt",
                            header=TRUE, sep='\t')
prot_stats = read.table("../prot_corrs/prot_stats.txt",
                        header=TRUE, sep='\t')
all_gene_tf_prot = sapply(all_genes, function(gene) {
  sapply(unique(prot_stats %>% filter(n_vals >= 20) %>% pull(tf) %>% as.character()), 
         function(tf) {
    paste(gene,tf,sep='_')
  }) })
fi_tsts_prot = do.call('rbind', 
                        lapply(unique(as.character(sig_hits$TISSUE_ABBRV)), function(tiss1) {
                          print(c(tiss1,'prot'))
                          this_table=table(all_gene_tf_prot %in% (sig_hits %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                           all_gene_tf_prot %in% (sig_hits_prot %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)))
                          this_fi = fisher.test(this_table)
                          print(this_fi)
                          data.frame(tiss1=c(tiss1,'prot'), tiss2=c('prot',tiss1),
                                     fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                        } ) )
fi_tsts_prot_mat = fi_tsts_prot %>%
  filter(tiss2=='prot') %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_prot_mat_mod = apply(fi_tsts_prot_mat, c(1,2), function(x) {
  max(min(x, 13), -4)
})

this_table=table(all_gene_tf_prot %in% (sig_hits %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                 all_gene_tf_prot %in% (sig_hits_prot %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)))
this_fi = fisher.test(this_table)
print(this_fi)
data.frame(tiss1=c('all_tiss'), tiss2=c('prot'),
           fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
.Machine$double.xmin

                       
fi_tst_cross_prot = table(all_gene_tf_prot %in% (sig_hits_cross %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)),
                          all_gene_tf_prot %in% (sig_hits_prot %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)) ) %>%
  fisher.test()
fi_tst_cross_prot

all_cross_tsts = rbind(fi_tsts_cross,fi_tsts_prot) %>%
  rbind(data.frame(tiss1=c('cross','cross','prot','prot'),
                   tiss2=c('cross','prot','cross','prot'),
                   fi_or=c(Inf,fi_tst_cross_prot$estimate,fi_tst_cross_prot$estimate,Inf),
                   fi_logor=log2(c(Inf,fi_tst_cross_prot$estimate,fi_tst_cross_prot$estimate,Inf)),
                   fi_p=c(0,fi_tst_cross_prot$p.value,fi_tst_cross_prot$p.value,0) ) )
all_cross_tsts_mat = all_cross_tsts %>%
  filter(tiss2 %in% c('cross','prot')) %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
all_cross_tsts_mat_mod = apply(all_cross_tsts_mat, c(1,2), function(x) {
  max(min(x, 13), -4)
})
heatmap3(all_cross_tsts_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2')


all_cross_tsts_winorder = all_cross_tsts %>%
  mutate(tiss1 = factor(tiss1, 
                        levels=c('cross','prot',
                                 'SPLEEN','LIVER','PNCREAS','PTTARY','NERVET','ARTTBL',
                                 'THYROID','ESPMSL','MSCLSK','ADPSBQ','LCL','LUNG','SKINS',
                                 'BRNCTXA','BRNCHA','BRNNCC'))) %>%
  arrange(tiss1)
all_cross_tsts_winorder_mat = all_cross_tsts_winorder %>%
  filter(tiss2 %in% c('cross','prot')) %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
all_cross_tsts_winorder_mat_mod = apply(all_cross_tsts_winorder_mat, c(1,2), function(x) {
  max(min(x, 13), -4)
})
heatmap3(all_cross_tsts_winorder_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2',
         Rowv = NA)





# #### testing out TF clustering...
# 
# 
# #all_genes = genes %>% pull(gene) %>% as.character()
# all_tiss = unique(as.character(sig_hits$TISSUE_ABBRV))
# all_gene_tiss = sapply(all_genes, function(gene_i) {
#   sapply(all_tiss, function(tiss_i) {
#     paste(gene_i,tiss_i,sep='_')
#   })
# })
# 
# # fi_tsts_fdr20_tfs = do.call('rbind', 
# #                         lapply(sort(unique(as.character(sig_hits$tf))), function(tf1) {
# #                           do.call('rbind',
# #                                   lapply(sort(unique(as.character(sig_hits$tf))), function(tf2) {
# #                                     print(c(tf1,tf2))
# #                                     this_table=table(all_gene_tiss %in% (sig_hits %>% filter(tf==tf1) %>% unite(gene_tiss,phenotype_id,TISSUE_ABBRV) %>% pull(gene_tiss)),
# #                                                      all_gene_tiss %in% (sig_hits %>% filter(tf==tf2) %>% unite(gene_tiss,phenotype_id,TISSUE_ABBRV) %>% pull(gene_tiss)))
# #                                     this_fi = fisher.test(this_table)
# #                                     print(this_fi)
# #                                     data.frame(tf1, tf2,
# #                                                fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
# #                                   }) )
# #                         } ) )
# # write.table(fi_tsts_fdr20_tfs, 'fi_tsts_fdr20_tfs.txt',
# #             col.names=TRUE, row.names=FALSE,
# #             sep='\t', quote = FALSE)
# fi_tsts_fdr20_tfs = read.table("fi_tsts_fdr20_tfs.txt",
#                                header=TRUE,sep='\t')
# 
# fi_tsts_fdr20_tfs_mat = fi_tsts_fdr20_tfs %>%
#   select(c(tf1,tf2,'fi_logor')) %>%
#   pivot_wider(names_from=tf2,
#               values_from=fi_logor) %>%
#   column_to_rownames("tf1") %>%
#   as.matrix()
# fi_tsts_fdr20_tfs_mat_mod = apply(fi_tsts_fdr20_tfs_mat, c(1,2), function(x) {
#   max(min(x, 12), -2)
# })
# heatmap3(fi_tsts_fdr20_tfs_mat_mod,
#          showColDendro = FALSE,
#          scale='none', balanceColor = TRUE,
#          method='ward.D2',
#          cexRow = .3, cexCol = .3)
# 
# 
# 
# 
# 
# 
