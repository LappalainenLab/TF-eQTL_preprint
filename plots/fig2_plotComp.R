#!/usr/bin/Rscript
##
##  fig2_plotComp.R
##
##  EDF 3/29/2021
##


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(heatmap3)

setwd("~/projects/MANUSCRIPT/")


tf_info = read.table("../TFi-eQTL/input_files/TFs.info.curated.txt",
                     header=TRUE, sep='\t')
tiss_info = read.table("../TFi-eQTL/input_files/tissue_info.final.txt",
                       header=TRUE, sep='\t')
genes = read.table("data_tables/variant_sets/caviar_chipseq_motif_variants.MAF05.genesort.txt",
                   header=TRUE, sep='\t') %>%
  group_by(gene) %>%
  summarize(num=n())


fi_tsts_fdr20 = read.table("data_tables/combine_corrs/fi_tsts_win_fdr20.txt",
                           header=TRUE,sep='\t')
fi_tsts_cross = read.table("data_tables/combine_corrs/fi_tsts_expr_prot.txt",
                           header=TRUE,sep='\t')


fi_tsts_fdr20_mat = fi_tsts_fdr20 %>%
  select(c(tiss1,tiss2,'fi_logor')) %>%
  pivot_wider(names_from=tiss2,
              values_from=fi_logor) %>%
  column_to_rownames("tiss1") %>%
  as.matrix()
fi_tsts_fdr20_mat_mod = apply(fi_tsts_fdr20_mat, c(1,2), function(x) {
  max(min(x, 13), 0)
})
pdf(file = "plots/fig2_win_comp.pdf",width=5,height=5)
heatmap3(fi_tsts_fdr20_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2')
dev.off()

all_cross_tsts_winorder = fi_tsts_cross %>%
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
  max(min(x, 13), -3)
})

pdf(file="plots/fig2_cross_comp.pdf",width=5,height=5)
heatmap3(all_cross_tsts_winorder_mat_mod,
         showColDendro = FALSE,
         scale='none', balanceColor = TRUE,
         method='ward.D2',
         Rowv = NA)
dev.off()



