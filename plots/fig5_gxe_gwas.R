#!/usr/bin/Rscript
##
##  fig5_gxe_gwas.R
##
##  EDF 7/9/21
##

setwd("~/projects/MANUSCRIPT/")

library(dplyr)
library(tidyr)
library(ggvenn)


rep_data = read.table("~/data/Findley2021/cASE_GxE_rep.tab",
                      header=TRUE,sep='\t')
sig_genes = read.table("~/projects/overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                       header=TRUE, sep='\t') %>%
  separate(phenotype_id, c('ensg'), sep='[.]',
           remove=FALSE)

genes_overlap_list = list(TF.eQTL = as.character(sig_genes$ensg),
                          GxE = as.character(rep_data$ensg))
ggvenn( genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GxE_gene.pdf',
       width=4, height=3)



moc_hits = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.vars_byset.txt",
                      header=TRUE,sep='\t') %>%
  unite(tiss_gene, corr, gene, remove=FALSE)
coloc_table = read.table("~/data/gtex/v8/GWAS_enloc/enloc_ENLOC_rcp_gt_0.5_with_gwas_pval.tsv",
                         header=TRUE,sep='\t') %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE)

genes_overlap_list = list(TF.eQTL = as.character(moc_hits$gene),
                          GWAS = as.character(coloc_table$gene_id))
ggvenn( genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GWAScoloc_gene.pdf',
       width=4, height=3)


all_eqtls = read.table('~/projects/examples/input_files/GTEx_Analysis_v8_eQTL.all_tissues.sig_egenes.txt.gz',
                       header=TRUE,sep='\t')
all_vars = read.table("~/projects/examples/input_files/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t')
all_eqtls_tested = all_eqtls %>%
  filter(gene_id %in% as.character(all_vars$gene)) %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE)
all_eqtls_tested_tiss = all_eqtls_tested %>%
  filter(tissue %in% as.character(moc_hits$corr))

tissue_genes_overlap_list = list(TF.eQTL = filter(moc_hits, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                                   pull(tiss_gene) %>% as.character(),
                                 GWAS = filter(coloc_table, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                                   pull(tiss_gene) %>% as.character())
ggvenn( tissue_genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('plots/fig5_GWAScoloc_tissgene.pdf',
       width=4, height=3)





