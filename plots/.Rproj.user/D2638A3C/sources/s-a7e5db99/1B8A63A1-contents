#!/usr/bin/Rscript
##
##  coloc_overlap.R
##
##  EDF 6/21/21
##

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(forcats)
library(gaston)
#library(VennDiagram)
library(ggvenn)

setwd("~/projects/examples/")

coloc_table = read.table("input_files/enloc_ENLOC_rcp_gt_0.5_with_gwas_pval.tsv",
                         header=TRUE,sep='\t') %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE)
moc_hits = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.vars_byset.txt",
                      header=TRUE,sep='\t') %>%
  unite(tiss_gene, corr, gene, remove=FALSE)
all_eqtls = read.table('input_files/GTEx_Analysis_v8_eQTL.all_tissues.sig_egenes.txt.gz',
                       header=TRUE,sep='\t')
all_vars = read.table("input_files/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t')


all_eqtls_tested = all_eqtls %>%
  filter(gene_id %in% as.character(all_vars$gene)) %>%
  unite(tiss_gene, tissue, gene_id, remove=FALSE)

## looking at tested genes, we see enrichment
table(unique(all_eqtls_tested$gene_id) %in% as.character(coloc_table$gene_id),
      unique(all_eqtls_tested$gene_id) %in% as.character(moc_hits$gene))
table(unique(all_eqtls_tested$gene_id) %in% as.character(coloc_table$gene_id),
             unique(all_eqtls_tested$gene_id) %in% as.character(moc_hits$gene)) %>% 
  fisher.test()
sum(unique(all_eqtls_tested$gene_id) %in% as.character(coloc_table$gene_id))/length(unique(all_eqtls_tested$gene_id))
sum(unique(all_eqtls_tested$gene_id) %in% as.character(coloc_table$gene_id) & unique(all_eqtls_tested$gene_id) %in% as.character(moc_hits$gene)) /
  sum(unique(all_eqtls_tested$gene_id) %in% as.character(moc_hits$gene))

genes_overlap_list = list(TF.eQTL = as.character(moc_hits$gene),
                          GWAS = as.character(coloc_table$gene_id))
ggvenn( genes_overlap_list,
                       fill_color=c('darkorchid4', 'gray') )
ggsave('GWAScoloc_gene_venn.pdf',
       width=4, height=3)

# table(all_eqtls_tested$gene_id %in% as.character(coloc_table$gene_id),
#       all_eqtls_tested$gene_id %in% as.character(moc_hits$gene))
# table(all_eqtls_tested$gene_id %in% as.character(coloc_table$gene_id),
#       all_eqtls_tested$gene_id %in% as.character(moc_hits$gene)) %>% 
#   fisher.test()
# sum(all_eqtls_tested$gene_id %in% as.character(coloc_table$gene_id))/nrow(all_eqtls_tested)
# sum(all_eqtls_tested$gene_id %in% as.character(coloc_table$gene_id) & all_eqtls_tested$gene_id %in% as.character(moc_hits$gene)) /
#   sum(all_eqtls_tested$gene_id %in% as.character(moc_hits$gene))

# looking at tested genes & tiss, and requiring moc tiss == coloc tiss
all_eqtls_tested_tiss = all_eqtls_tested %>%
  filter(tissue %in% as.character(moc_hits$corr))
table(unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(coloc_table$tiss_gene),
      unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(moc_hits$tiss_gene))
table(unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(coloc_table$tiss_gene),
      unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(moc_hits$tiss_gene)) %>%
  fisher.test()
sum(unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(coloc_table$tiss_gene))/length(unique(all_eqtls_tested_tiss$tiss_gene))
sum(unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(coloc_table$tiss_gene) & unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(moc_hits$tiss_gene)) /
  sum(unique(all_eqtls_tested_tiss$tiss_gene) %in% as.character(moc_hits$tiss_gene))


tissue_genes_overlap_list = list(TF.eQTL = filter(moc_hits, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                                   pull(tiss_gene) %>% as.character(),
                          GWAS = filter(coloc_table, tiss_gene %in% all_eqtls_tested_tiss$tiss_gene) %>% 
                            pull(tiss_gene) %>% as.character())
ggvenn( tissue_genes_overlap_list,
        fill_color=c('darkorchid4', 'gray') )
ggsave('GWAScoloc_tissgene_venn.pdf',
       width=4, height=3)

all_eqtls_tested_tiss2 = all_eqtls_tested %>%
  filter(tissue %in% as.character(moc_hits$corr)) %>%
  group_by(gene_id) %>%
  summarize(n_rows = n(),
            coloc_overlap = sum(tiss_gene %in% coloc_table$tiss_gene),
            moc_overlap = sum(tiss_gene %in% moc_hits$tiss_gene))

table(all_eqtls_tested_tiss2$coloc_overlap > 0,
      all_eqtls_tested_tiss2$moc_overlap > 0)
table(all_eqtls_tested_tiss2$coloc_overlap > 0,
      all_eqtls_tested_tiss2$moc_overlap > 0) %>%
  fisher.test()
sum(all_eqtls_tested_tiss2$coloc_overlap > 0)/nrow(all_eqtls_tested_tiss2)
sum(all_eqtls_tested_tiss2$coloc_overlap > 0 & all_eqtls_tested_tiss2$moc_overlap > 0) / 
  sum(all_eqtls_tested_tiss2$moc_overlap > 0)



all_eqtls_tested_tiss3 = all_eqtls_tested %>%
  filter(tissue %in% as.character(moc_hits$corr)) %>%
  mutate(coloc_overlap_tg = tiss_gene %in% coloc_table$tiss_gene,
         moc_overlap_tg = tiss_gene %in% moc_hits$tiss_gene) %>%
  group_by(gene_id) %>%
  summarize(n_rows = n(),
            coloc_overlap = sum(coloc_overlap_tg),
            moc_overlap = sum(moc_overlap_tg),
            both_overlap = sum(coloc_overlap_tg > 0 & moc_overlap_tg > 0))
sum(all_eqtls_tested_tiss3$coloc_overlap > 0)/nrow(all_eqtls_tested_tiss3)
sum(all_eqtls_tested_tiss3$both_overlap > 0) / 
  sum(all_eqtls_tested_tiss3$moc_overlap > 0)





gwas_info = read.table("input_files/gwas/GWAS_metadata_TableS11_PMID32913098.abbrv.csv",
                       header=TRUE, sep='\t')

coloc_moc_hits = all_eqtls_tested %>%
  merge(coloc_table,
        by='tiss_gene') %>%
  merge(moc_hits,
        by='tiss_gene') %>%
  merge(gwas_info,
        by.x='phenotype',
        by.y='Tag') %>%
  select(phenotype, Category, gene, gene_name, tf, tissue.x, rcp, lead_snp, variant, variant_id) %>%
  rename(tissue=tissue.x, enloc_rcp = rcp, lead_coloc_snp = lead_snp, top_tf_var = variant, top_eqtl_var = variant_id) %>%
  arrange(phenotype, gene, tf, tissue)


bed_matrix_gts = read.bed.matrix("/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/plink_genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01")
coloc_moc_hits$ld = apply(coloc_moc_hits, 1, function(rowi) {
  tf_var = rowi[9]
  coloc_var = rowi[8]
  LD(bed_matrix_gts,
     which(bed_matrix_gts@snps$id == tf_var),
     which(bed_matrix_gts@snps$id == coloc_var))
})

coloc_moc_hits_filt = coloc_moc_hits %>%
  filter(ld > 0.4) %>%
  mutate(Category = ifelse(Category == 'Allergy', 'Immune', 
                           as.character(Category)))

write.table(coloc_moc_hits_filt,
            "coloc_moc_hits.txt",
            col.names=TRUE,sep='\t',
            row.names=FALSE,quote=FALSE)


nrow(coloc_moc_hits_filt)
coloc_moc_hits_filt %>%
  unite(gene_gwas_tiss, gene_name, phenotype, tissue) %>%
  pull(gene_gwas_tiss) %>%
  unique() %>%
  length()
coloc_moc_hits_filt %>%
  unite(gene_gwas, gene_name, phenotype) %>%
  pull(gene_gwas) %>%
  unique() %>%
  length()


table(coloc_moc_hits_filt %>%  
        group_by(phenotype, tf, gene) %>%
        summarize(n=n()) %>%
        unite(ptf, phenotype, tf) %>% 
        select(ptf)) %>% sort()

coloc_moc_hits_filt %>%  
  group_by(Category, tf, gene) %>%
  summarize(n=n()) %>%
  group_by(Category, tf) %>%
  summarize(n=n())

coloc_moc_hits_filt %>%  
  group_by(Category, tf, gene) %>%
  summarize(n=n()) %>%
  group_by(Category, tf) %>%
  summarize(n=n()) %>%
  group_by(Category) %>%
  summarize(max_n = max(n),
            max_tfs = tf[which(n==max_n)]) %>%
  View()


coloc_moc_hits_filt %>%  
  group_by(Category, tf) %>%
  summarize(n_hits=n(),
            n_gene = length(unique(gene))) %>%
  filter(n_gene>2) %>% 
  View()

coloc_moc_hits_filt %>%
  group_by(Category, tf) %>%
  summarize(n_hits=n(),
            n_gene = length(unique(gene))) %>%
  filter(n_gene>1,
         ifelse(Category=='Anthropometric', n_gene>4,
                ifelse(Category=='Blood',n_gene>3,
                       n_gene>1))) %>%
  ggplot(aes(Category, n_gene)) +
  geom_col(aes(group=tf,
               fill=tf),
           position=position_dodge(width=.9)) +
  geom_text(aes(group=tf,
                 label=tf),
             position=position_dodge(width=.9),
            angle=90, hjust=0) +
  ylim(0,7) +
  scale_fill_viridis_d() +
  theme_classic()

coloc_moc_hits_filt %>%
  filter(Category == 'Immune',
         tf=='CTCF')


coloc_moc_hits_filt %>%
  filter(Category == 'Anthropometric',
         tf=='RXRA')

coloc_moc_hits_filt %>%
  filter(Category == 'Blood',
         tf=='MAX')


coloc_moc_hits_filt %>%
  filter(Category == 'Psychiatric_neurologic',
         tf=='STAT3')



# 
# 
# table(all_eqtls_tested_tiss$tiss_gene %in% as.character(coloc_table$tiss_gene),
#       all_eqtls_tested_tiss$tiss_gene %in% as.character(moc_hits$tiss_gene))
# table(all_eqtls_tested_tiss$tiss_gene %in% as.character(coloc_table$tiss_gene),
#       all_eqtls_tested_tiss$tiss_gene %in% as.character(moc_hits$tiss_gene)) %>% 
#   fisher.test()
# sum(all_eqtls_tested_tiss$tiss_gene %in% as.character(coloc_table$tiss_gene))/nrow(all_eqtls_tested_tiss)
# sum(all_eqtls_tested_tiss$tiss_gene %in% as.character(coloc_table$tiss_gene) & all_eqtls_tested_tiss$tiss_gene %in% as.character(moc_hits$tiss_gene)) / 
#   sum(all_eqtls_tested_tiss$tiss_gene %in% as.character(moc_hits$tiss_gene))
# 



coloc_moc_hits2 = all_eqtls_tested %>%
  merge(coloc_table,
        by=c('tiss_gene','gene_id','tissue')) %>%
  merge(moc_hits,
        by.x='gene_id',
        by.y='gene') %>%
  merge(gwas_info,
        by.x='phenotype',
        by.y='Tag') %>%
  select(phenotype, Category, gene_id, gene_name, tf, rcp, lead_snp, variant, variant_id, tissue, corr) %>%
  rename(gene=gene_id, enloc_rcp = rcp, 
         lead_coloc_var = lead_snp, top_tf_var = variant, top_eqtl_var = variant_id,
         coloc_tiss = tissue, tfeqtl_tiss = corr) %>%
  arrange(phenotype, gene, tf)
coloc_moc_hits2$ld = apply(coloc_moc_hits2, 1, function(rowi) {
  tf_var = rowi[8]
  coloc_var = rowi[7]
  LD(bed_matrix_gts,
     which(bed_matrix_gts@snps$id == tf_var),
     which(bed_matrix_gts@snps$id == coloc_var))
})

write.table(coloc_moc_hits2,
            "coloc_moc_hits2.txt",
            col.names=TRUE,sep='\t',
            row.names=FALSE,quote=FALSE)

coloc_moc_hits2_filt = coloc_moc_hits2 %>%
  filter(ld > 0.4) %>%
  mutate(Category = ifelse(Category == 'Allergy', 'Immune', 
                           as.character(Category)))


coloc_moc_hits2_filt %>%
  group_by(Category, tf) %>%
  summarize(n_hits=n(),
            n_gene = length(unique(gene))) %>%
  filter(n_gene>2,
         ifelse(Category=='Anthropometric', n_gene>9,
                ifelse(Category=='Blood',n_gene>8,
                       ifelse(Category %in% c('Cardiometabolic','Psychiatric_neurologic'),n_gene>3,
                       n_gene>1)))) %>%
  ggplot(aes(Category, n_gene)) +
  geom_col(aes(group=tf,
               fill=tf),
           position=position_dodge(width=.9)) +
  geom_text(aes(group=tf,
                label=tf),
            position=position_dodge(width=.9),
            angle=90, hjust=0) +
  ylim(0,12) +
  scale_fill_viridis_d() +
  theme_classic()

coloc_moc_hits2_filt %>%
  group_by(phenotype, tf) %>%
  summarize(n_hits=n(),
            n_gene = length(unique(gene))) %>%
  filter(n_gene>3) %>%
  ggplot(aes(phenotype, n_gene)) +
  geom_col(aes(group=tf,
               fill=tf),
           position=position_dodge(width=.9)) +
  geom_text(aes(group=tf,
                label=tf),
            position=position_dodge(width=.9),
            angle=90, hjust=0) +
  ylim(0,6) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))



coloc_moc_hits2_filt %>%
  filter(Category == 'Anthropometric',
         tf=='CTCF') %>%
  group_by(phenotype,gene_name) %>%
  summarize(n=n())
