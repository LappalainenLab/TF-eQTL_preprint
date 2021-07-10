#!/usr/bin/Rscript
##
##  analyzeFishersIRF1ASE_allgenes_timecomb.R
##
##  EDF 4/15/21
##


library(dplyr)
library(ggplot2)
library(tidyr)
library(metap)
library(UpSetR)

setwd("~/projects/IRF1/")



mpileups_comb_fi_test = read.table("results/Fi_tst_combined_byvar_timecomb.allgenes.txt",
                                         header=TRUE, sep='\t')
# mpileups_comb_fi_test_bygene = read.table("results/Fi_tst_combined_bygene.allgenes.txt",
#                                          header=TRUE, sep='\t')

moc_tf_eqtls = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                          header=TRUE,sep='\t')

table(unique(mpileups_comb_fi_test$gene) %in% unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)),
      unique(mpileups_comb_fi_test$gene) %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
table(unique(mpileups_comb_fi_test$gene) %in% unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)),
      unique(mpileups_comb_fi_test$gene) %in% unique(moc_tf_eqtls %>% filter(startsWith(as.character(tf),'IRF')) %>% pull(phenotype_id)))
# 
# table(unique(mpileups_comb_fi_test_bygene$gene) %in% unique(mpileups_comb_fi_test_bygene %>% filter(fi_p < 0.05) %>% pull(gene)),
#       unique(mpileups_comb_fi_test_bygene$gene) %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
# table(unique(mpileups_comb_fi_test_bygene$gene) %in% unique(mpileups_comb_fi_test_allps_bygene %>% filter(fi_p < 0.05) %>% pull(gene)),
#       unique(mpileups_comb_fi_test_allps_bygene$gene) %in% unique(moc_tf_eqtls %>% filter(startsWith(as.character(tf),'IRF')) %>% pull(phenotype_id)))






eqtl_snps = read.table("eqtls/gene_var_gt.multi_or_cross.hets.txt",
                       header=FALSE, sep='\t') %>%
  mutate(chr_pos = paste(V3,V4,sep='_'))
IRF1_top_snps = read.table("input_files/TFi_hits_topvars.IRF1.bed",
                      header=FALSE,sep='\t') %>%
  mutate(chr_pos=paste(V1,V3,sep = '_'))

top_snps_IRF1_het_in_HEK = eqtl_snps %>% filter(chr_pos %in% IRF1_top_snps$chr_pos)

table(unique(mpileups_comb_fi_test$gene) %in% unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)),
      unique(mpileups_comb_fi_test$gene) %in% unique(top_snps_IRF1_het_in_HEK$V1)) #%>% fisher.test()

# table(unique(mpileups_comb_fi_test_allps_bygene$gene) %in% unique(mpileups_comb_fi_test_bygene %>% filter(fi_p < 0.05) %>% pull(gene)),
#       unique(mpileups_comb_fi_test_bygene$gene) %in% unique(top_snps_IRF1_het_in_HEK$V1)) #%>% fisher.test()


mpileups_comb_fi_test %>% filter(gene %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
mpileups_comb_fi_test %>% filter(gene %in% unique(top_snps_IRF1_het_in_HEK$V1))

# mpileups_comb_fi_test_bygene %>% filter(gene %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
# mpileups_comb_fi_test_bygene %>% filter(gene %in% unique(top_snps_IRF1_het_in_HEK$V1))




qqplot(-log10(mpileups_comb_fi_test$fi_p),
        -log10(mpileups_comb_fi_test$chi_p))
abline(a=0,b=1)

# qqplot(-log10(mpileups_comb_fi_test$fi_p),
#        -log10(mpileups_comb_fi_test_bygene %>% 
#                 filter(gene %in% unique(top_snps_IRF1_het_in_HEK$V1)) %>%
#                 pull(fi_p)))
# abline(a=0,b=1)



# 
# top_eqtl_snps = merge(top_snps,eqtl_snps,by='chr_pos')
# length(unique(top_eqtl_snps$V1.y))
# 
# temp = read.table("coding_snps/gene_exon_var_gt.multi_or_cross.hets.txt",header=TRUE,sep='\t')
# table(unique(temp$gene) %in% as.character(unique(top_eqtl_snps$V1.y)))
# table(as.character(unique(top_eqtl_snps$V1.y)) %in% unique(temp$gene))
# 
# 
# tiss_hits_fdr20 = read.table("../TFi-eQTL/tensorqtl/sig_assoc.fdr20.partial.txt",
#                              header=TRUE, sep='\t')






## REGULON ENRICHMENT

irf1_reg = read.table("input_files/merged_interactions.txt",
                      header=TRUE,sep='\t') %>%
  filter(TF=='IRF1')

irf_reg = read.table("input_files/merged_interactions.txt",
                     header=TRUE,sep='\t') %>%
  filter(startsWith(as.character(TF),'IRF'))

stat_reg = read.table("input_files/merged_interactions.txt",
                      header=TRUE,sep='\t') %>%
  filter(startsWith(as.character(TF),'STAT'))

nfkb_reg = read.table("input_files/merged_interactions.txt",
                      header=TRUE,sep='\t') %>%
  filter(startsWith(as.character(TF),'NFKB'))



table(unique(mpileups_comb_fi_test$description) %in% unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(description)),
      unique(mpileups_comb_fi_test$description) %in% unique(as.character(irf1_reg$target)))
table(unique(mpileups_comb_fi_test$description) %in% unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(description)),
      unique(mpileups_comb_fi_test$description) %in% unique(as.character(irf_reg$target)))
# 
# table(unique(mpileups_comb_fi_test_bygene$description) %in% unique(mpileups_comb_fi_test_bygene %>% filter(fi_p < 0.05) %>% pull(description)),
#       unique(mpileups_comb_fi_test_bygene$description) %in% unique(as.character(irf1_reg$target)))
# table(unique(mpileups_comb_fi_test_bygene$description) %in% unique(mpileups_comb_fi_test_bygene %>% filter(fi_p < 0.05) %>% pull(description)),
#       unique(mpileups_comb_fi_test_bygene$description) %in% unique(as.character(irf_reg$target)))






irf1_genes_cross = read.table("../overlap/replication/sig_assoc.fdr05.cross_hits.txt",
                              header=TRUE,sep='\t') %>%
  filter(tf=='IRF1')
irf1_genes_tiss = read.table("../overlap/replication/sig_assoc.fdr20.any_tiss_hits.txt",
                             header=TRUE,sep='\t') %>%
  filter(tf=='IRF1')

all_tested_genes = unique(as.character(mpileups_comb_fi_test$gene))

upset(fromList(list(
  cross = irf1_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
    pull(phenotype_id) %>% as.character(),
  tiss = irf1_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>%
    pull(phenotype_id) %>% as.character(),
  HEK = mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene) %>% as.character())))

fisher.test(table(all_tested_genes %in% as.character((irf1_genes_cross$phenotype_id)),
                  all_tested_genes %in% as.character((irf1_genes_tiss$phenotype_id))))
fisher.test(table(all_tested_genes %in% (mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)),
                  all_tested_genes %in% as.character((irf1_genes_tiss$phenotype_id))))
fisher.test(table(all_tested_genes %in% (mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)),
                  all_tested_genes %in% as.character((irf1_genes_cross$phenotype_id))))


mpileups_comb_fi_test %>%
  pull(fi_p) %>%
  hist()

mpileups_comb_fi_test %>%
  filter(gene %in% (irf1_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
                      pull(phenotype_id) %>% as.character())) %>%
  pull(fi_p) %>%
  hist()

mpileups_comb_fi_test %>%
  filter(gene %in% (irf1_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>% 
                      pull(phenotype_id) %>% as.character())) %>%
  pull(fi_p) %>%
  hist()

qqplot(mpileups_comb_fi_test %>%
         pull(fi_p) %>%
         log10() * -1,
       mpileups_comb_fi_test %>%
         filter(gene %in% (irf1_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
                             pull(phenotype_id) %>% as.character())) %>%
         pull(fi_p) %>%
         log10() * -1)
abline(a=0,b=1)

qqplot(mpileups_comb_fi_test %>%
         pull(fi_p) %>%
         log10() * -1,
       mpileups_comb_fi_test %>%
         filter(gene %in% (irf1_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>% 
                             pull(phenotype_id) %>% as.character())) %>%
         pull(fi_p) %>%
         log10() * -1)
abline(a=0,b=1)

qqplot(mpileups_comb_fi_test %>%
         pull(fi_p) %>%
         log10() * -1,
       mpileups_comb_fi_test %>%
         filter(gene %in% (irf1_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>% 
                             pull(phenotype_id) %>% as.character()) |
                  gene %in% (irf1_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
                               pull(phenotype_id) %>% as.character())) %>%
         pull(fi_p) %>%
         log10() * -1)
abline(a=0,b=1)



irf1_ps_cross = read.table("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.IRF1.top.adjp.txt",
                           header=TRUE, sep='\t')
qqplot(irf1_ps_cross %>%
         filter(gene %in% all_tested_genes) %>%
         pull(sp_p.bf_meff_gene.min) %>% log10() * -1,
       irf1_ps_cross %>%
         filter(gene %in% (mpileups_comb_fi_test %>%
                             filter(fi_p < 0.05) %>%
                             pull(gene) %>% as.character())) %>%
         pull(sp_p.bf_meff_gene.min) %>% log10() * -1)
abline(a=0,b=1)


irf1_tiss_files = unlist(lapply(list.files("../TFi-eQTL/tensorqtl/","*",full.names=TRUE), function(folderi) {
  print(paste0(folderi,"/IRF1/"))
  list.files(paste0(folderi,"/IRF1/"),
             "*.IRF1.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz",
                              full.names = TRUE)
}))
irf1_ps_tiss = do.call('rbind', lapply(irf1_tiss_files, function(filei) {
  tissi = strsplit(filei,"/")[[1]][5]
  read.table(filei, header=TRUE, sep='\t') %>%
    mutate(tiss=tissi)
} ) ) %>%
  filter(! (tiss %in% 
              c('Whole_Blood','Colon_Transverse','Cells_Cultured_fibroblasts','Stomach','Testis')))
qqplot(irf1_ps_tiss %>%
         filter(phenotype_id %in% all_tested_genes) %>%
         pull(pval_emt) %>% log10() * -1,
       irf1_ps_tiss %>%
         filter(phenotype_id %in% (mpileups_comb_fi_test %>%
                             filter(fi_p < 0.05) %>%
                             pull(gene) %>% as.character())) %>%
         pull(pval_emt) %>% log10() * -1)
abline(a=0,b=1)






irfall_genes_cross = read.table("../overlap/replication/sig_assoc.fdr05.cross_hits.txt",
                              header=TRUE,sep='\t') %>%
  filter(startsWith(as.character(tf),'IRF'))
irfall_genes_tiss = read.table("../overlap/replication/sig_assoc.fdr20.any_tiss_hits.txt",
                             header=TRUE,sep='\t') %>%
  filter(startsWith(as.character(tf),'IRF'))

upset(fromList(list(
  cross = irfall_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
    pull(phenotype_id) %>% as.character() %>% unique(),
  tiss = irfall_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>%
    pull(phenotype_id) %>% as.character() %>% unique(),
  HEK = mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene) %>% as.character() %>% unique())))



## Look at top examps not from MOC
mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% arrange(fi_p)

mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% filter(description %in% irf1_reg$target)

irf1_overlap = read.table("../crosstiss_tf_corrs/overlaps/IRF1.txt.gz",
                          header=TRUE,sep='\t')
mpileups_comb_fi_test %>%
  filter(fi_p < 0.05) %>%
  filter(description %in% 
           (irf1_overlap %>%
              filter(tf_chip+tf_motif > 0) %>%
              pull(gene_name)))

irf1_cross_hits = read.table("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.IRF1.top.adjp.txt",
                             header=TRUE,sep='\t') %>%
  filter(sp_p.bf_meff_gene.min.bh_tf < 0.05)
irf1_any_tiss_hits = do.call('rbind', lapply(list.files("../overlap/replication/indiv_tiss/",
                                                        "*fdr20*",
                                                        full.names = TRUE), 
                                             function(filei) {
                                               read.table(filei,
                                                          header=TRUE,sep='\t') %>%
                                                 filter(tf=='IRF1')
                                             } ))
mpileups_comb_fi_test %>%
  filter(fi_p < 0.05) %>%
  filter(description %in% 
           (irf1_overlap %>%
              filter(tf_chip+tf_motif > 0) %>%
              pull(gene_name)),
         gene %in% c(as.character(irf1_cross_hits$gene), 
                     as.character(irf1_any_tiss_hits$phenotype_id))) %>%
  arrange(fi_p)




table(unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_cross_hits$gene),
      unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_any_tiss_hits$phenotype_id))
table(unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_genes_cross$phenotype_id),
      unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_genes_tiss$phenotype_id))


upset(fromList(list(
  cross = irf1_genes_cross %>% filter(phenotype_id %in% all_tested_genes) %>% 
    pull(phenotype_id) %>% as.character(),
  tiss = irf1_genes_tiss %>% filter(phenotype_id %in% all_tested_genes) %>%
    pull(phenotype_id) %>% as.character(),
  HEK = mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene) %>% as.character())))



HEK_gts = read.table("input_files/293_CG.liftoverhg38.sorted.cleaned.vcf.gz",
                     header=FALSE,sep='\t')
HEK_gts_het = HEK_gts %>%
  filter(V10=='0/1') %>%
  unite(var,c(V1,V2,V4,V5),sep='_') %>%
  mutate(var_id = paste(var,'b38',sep='_'))


irf1_cross_gthet = merge(irf1_cross_hits, HEK_gts_het,
                      by.x = 'var.min', by.y='var_id')
irf1_win_gthet = merge(irf1_any_tiss_hits, HEK_gts_het,
                       by.x = 'variant_id', by.y = 'var_id')

table(unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_cross_gthet$gene),
      unique(mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene)) %in% 
        as.character(irf1_win_gthet$phenotype_id))


upset(fromList(list(
  cross = irf1_cross_gthet %>% filter(gene %in% all_tested_genes) %>% 
    pull(gene) %>% as.character(),
  tiss = irf1_win_gthet %>% filter(phenotype_id %in% all_tested_genes) %>%
    pull(phenotype_id) %>% as.character(),
  HEK = mpileups_comb_fi_test %>% filter(fi_p < 0.05) %>% pull(gene) %>% as.character())))






