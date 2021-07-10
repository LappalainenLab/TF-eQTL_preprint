#!/usr/bin/Rscript
##
##  plotTFeQTLs_compfig.R
##
##  EDF 1/16/20
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


sig_hits_tiss = read.table("../TFi-eQTL/eqtl_sig/all.fdr20.sig_eqtls.txt",
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
                        lapply(unique(as.character(sig_hits_tiss$TISSUE_ABBRV)), function(tiss1) {
                          do.call('rbind',
                                  lapply(unique(as.character(sig_hits_tiss$TISSUE_ABBRV)), function(tiss2) {
                                    print(c(tiss1,tiss2))
                                    this_table=table(all_gene_tf %in% (sig_hits_tiss %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                                     all_gene_tf %in% (sig_hits_tiss %>% filter(TISSUE_ABBRV==tiss2) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)))
                                    this_fi = fisher.test(this_table)
                                    print(this_fi)
                                    data.frame(tiss1, tiss2,
                                               fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                                  }) )
                        } ) )

write.table(fi_tsts_fdr20,
            "../TFi-eQTL/comp/fi_tsts_win_fdr20_new.txt",
            col.names=TRUE, sep='\t',
            quote=FALSE, row.names=FALSE)

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



sig_hits_cross = read.table("../crosstiss_tf_corrs/correlations/cross_tiss_tf_expr_corrs.med0.curated_set.fdr05.txt",
                            header=TRUE, sep='\t')
fi_tsts_cross = do.call('rbind', 
                        lapply(unique(as.character(sig_hits_tiss$TISSUE_ABBRV)), function(tiss1) {
                          print(c(tiss1,'cross'))
                          this_table=table(all_gene_tf %in% (sig_hits_tiss %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                           all_gene_tf %in% (sig_hits_cross %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)))
                          this_fi = fisher.test(this_table)
                          print(this_fi)
                          data.frame(tiss1=c(tiss1), tiss2=c('cross'),
                                     fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                        } ) )
write.table(fi_tsts_cross,
            "../TFi-eQTL/comp/fi_tsts_cross_new.txt",
            col.names=TRUE, sep='\t',
            quote=FALSE, row.names=FALSE)

fi_tsts_cross_mat = fi_tsts_cross %>%
  filter(tiss2=='cross') %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_cross_mat_mod = apply(fi_tsts_cross_mat, c(1,2), function(x) {
  max(min(x, 14), -4)
})


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
                        lapply(unique(as.character(sig_hits_tiss$TISSUE_ABBRV)), function(tiss1) {
                          print(c(tiss1,'prot'))
                          this_table=table(all_gene_tf_prot %in% (sig_hits_tiss %>% filter(TISSUE_ABBRV==tiss1) %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                                           all_gene_tf_prot %in% (sig_hits_prot %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)))
                          this_fi = fisher.test(this_table)
                          print(this_fi)
                          data.frame(tiss1=c(tiss1,'prot'), tiss2=c('prot',tiss1),
                                     fi_or = this_fi$estimate, fi_logor = log2(this_fi$estimate), fi_p = this_fi$p.value)
                        } ) )
write.table(fi_tsts_prot,
            "../TFi-eQTL/comp/fi_tsts_prot_new.txt",
            col.names=TRUE, sep='\t',
            quote=FALSE, row.names=FALSE)

fi_tsts_prot_mat = fi_tsts_prot %>%
  filter(tiss2=='prot') %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_prot_mat_mod = apply(fi_tsts_prot_mat, c(1,2), function(x) {
  max(min(x, 14), -4)
})

                       
fi_tst_cross_prot = table(all_gene_tf_prot %in% (sig_hits_cross %>% unite(gene_tf,gene,tf) %>% pull(gene_tf)),
                          all_gene_tf_prot %in% (sig_hits_prot %>% unite(gene_tf,phenotype_id,tf) %>% pull(gene_tf)) ) %>%
  fisher.test()


all_cross_tsts = rbind(fi_tsts_cross,fi_tsts_prot) %>%
  rbind(data.frame(tiss1=c('cross','cross','prot','prot'),
                   tiss2=c('cross','prot','cross','prot'),
                   fi_or=c(Inf,fi_tst_cross_prot$estimate,fi_tst_cross_prot$estimate,Inf),
                   fi_logor=log2(c(Inf,fi_tst_cross_prot$estimate,fi_tst_cross_prot$estimate,Inf)),
                   fi_p=c(0,fi_tst_cross_prot$p.value,fi_tst_cross_prot$p.value,0) ) )
write.table(all_cross_tsts,
            "../TFi-eQTL/comp/fi_tsts_expr_prot_new.txt",
            col.names=TRUE, sep='\t',
            quote=FALSE, row.names=FALSE)
all_cross_tsts_mat = all_cross_tsts %>%
  filter(tiss2 %in% c('cross','prot')) %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
all_cross_tsts_mat_mod = apply(all_cross_tsts_mat, c(1,2), function(x) {
  max(min(x, 14), -4)
})
heatmap3(all_cross_tsts_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2')


all_cross_tsts_winorder = all_cross_tsts %>%
  mutate(tiss1 = factor(tiss1, 
                        levels=c('cross','prot',
                                 'SPLEEN','LIVER','PNCREAS','PTTARY','NERVET','ARTTBL',
                                 'THYROID','ESPMSL','MSCLSK',
                                 'ADPSBQ','LCL','LUNG','SKINS',
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
  max(min(x, 14), -5)
})
heatmap3(all_cross_tsts_winorder_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2',
         Rowv = NA)





# any_hits = sig_hits_tiss %>%
#   group_by(tf,phenotype_id) %>%
#   summarize(num_tiss = n())
# nrow(any_hits)
# write.table(any_hits, "replication/sig_assoc.fdr20.any_tiss_hits.txt",
#             col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

multi_tiss_hits = sig_hits_tiss %>%
  group_by(tf,gene) %>%
  summarize(variants = paste(as.character(var),collapse=","),
            datasets = paste(as.character(tiss),collapse=','),
            num_tiss = n()) %>%
  filter(num_tiss > 1)
nrow(multi_tiss_hits)
table(multi_tiss_hits$num_tiss)
write.table(multi_tiss_hits, "../TFi-eQTL/comp/multi_within_tiss_corrs_fdr20.txt",
            col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)


cross_win_hits = rbind(sig_hits_tiss %>% select(tf,gene,var,tiss), 
                        sig_hits_cross %>% select(tf,gene, var) %>% mutate(tiss='cross')) %>%
  group_by(tf,gene) %>%
  summarize(variants = paste(as.character(var),collapse=","),
            datasets = paste(as.character(tiss),collapse=','),
            num_data = n(),
            cross = ifelse('cross' %in% tiss, TRUE,FALSE)) %>%
  filter(num_data > 1, cross) %>%
  select(-cross)
nrow(cross_win_hits)
write.table(cross_win_hits, "../TFi-eQTL/comp/within_cross_corrs.txt",
            col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)


multi_data_hits = rbind(sig_hits_tiss %>% select(tf,gene,var,tiss), 
                  sig_hits_cross %>% select(tf,gene, var) %>% mutate(tiss='cross')) %>%
  group_by(tf,gene) %>%
  summarize(variants = paste(as.character(var),collapse=","),
            datasets = paste(as.character(tiss),collapse=','),
            num_data = n()) %>%
  filter(num_data > 1)
nrow(multi_data_hits)
write.table(multi_data_hits, "../TFi-eQTL/comp/multi_or_with_cross_corrs.txt",
            col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)


# 
# multi_cross_hits = both_hits %>%
#   filter(num_tiss > 1)
# 
# nrow(multi_cross_hits)
# write.table(multi_cross_hits, "replication/sig_assoc.fdr20.multi_cross_hits.txt",
#             col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
# 
# multi_or_cross_hits = merge(both_hits, multi_tiss_hits,
#                             all = TRUE)
# nrow(multi_or_cross_hits)
# write.table(multi_or_cross_hits, "replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
#             col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)




