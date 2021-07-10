#!/usr/bin/Rscript
##
##  plotExamp1.R
##
##  EDF 6/9/2021
##

library(dplyr)
library(ggplot2)
library(gaston)
library(seqminer)
library(tidyr)
library(tibble)
library(ggsci)
library(scales)
library(ggrepel)

setwd("~/projects/examples/")

## TF: IKZF1
tfi='IKZF1'
## Gene: APBB1IP
descriptioni="APBB1IP"
## Trait: Blood cell counts

overlap_var = 'chr10_26438309_T_A_b38'


gene_info = read.table("input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                        header=TRUE,sep='\t')
genei=gene_info %>%
  filter(description==descriptioni) %>%
  pull(gene) %>%
  as.character()

tiss_info = read.table("input_files/tissue_translation_colors_v8.txt",
                       header=TRUE,sep='\t')

## tf-eqtls
cross_tiss_stats = read.table(paste0("input_files/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.",tfi,".top.adjp.txt"),
                              header=TRUE,sep='\t') %>%
  filter(gene == genei)

sig_tiss_stats = do.call('rbind', lapply(list.files("input_files/indiv_tiss/",
                                                    "*fdr20*",
                                                    full.names = TRUE), 
                                         function(filei) {
                                           read.table(filei,
                                                      header=TRUE,sep='\t') %>%
                                             filter(tf==tfi,
                                                    phenotype_id == genei)
                                         } ))
sig_tiss = sig_tiss_stats %>% pull(tiss) %>% as.character()
sig_tiss_short = tiss_info %>%
  filter(TISSUE_NAME %in% sig_tiss) %>%
  pull(TISSUE_ABBRV)

cross_tiss_var = cross_tiss_stats %>% pull(var.min) %>% as.character()
sig_tiss_var = sig_tiss_stats %>% pull(variant_id) %>% as.character()
print(paste("Cross tiss top var is",cross_tiss_var))
print(cross_tiss_stats$sp_p.bf_meff_gene.min.bh_tf)
print(paste(sig_tiss_stats %>% pull(tiss),"top var is",sig_tiss_var))


## overlap var
chri=strsplit(overlap_var,"_")[[1]][1]
posi=as.numeric(strsplit(overlap_var,"_")[[1]][2])


bed_matrix_gts = read.bed.matrix("/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/plink_genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01")
# lapply(sig_tiss_var, function(vari) {
#   LD(bed_matrix_gts,
#    which(bed_matrix_gts@snps$id == cross_tiss_var), 
#    which(bed_matrix_gts@snps$id == vari))
# })
lapply(c(cross_tiss_var, sig_tiss_var), function(vari) {
  LD(bed_matrix_gts,
   which(bed_matrix_gts@snps$id == overlap_var),
   which(bed_matrix_gts@snps$id == vari))
})

all_ld = LD(bed_matrix_gts,
            which(bed_matrix_gts@snps$id == overlap_var),
            c(which(bed_matrix_gts@snps$id == overlap_var)-1000,
              which(bed_matrix_gts@snps$id == overlap_var)+1000)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('var1') %>%
  rename(r2=chr10_26438309_T_A_b38) %>%
  mutate(r2_cat = factor(ifelse(var1 == overlap_var, 'lead',
                                ifelse(r2>.8,'.8<r2<=1',
                                       ifelse(r2>.6,'.6<r2<=.8',
                                              ifelse(r2>.4,'.4<r2<=.6',
                                                     ifelse(r2>.2,'.2<r2<=.4','r2<=.2'))))),
                         levels=c('.8<r2<=1','.6<r2<=.8','.4<r2<=.6','.2<r2<=.4','r2<=.2','lead')))


## coloc info
enloc_hits = read.table("input_files/enloc_ENLOC_rcp_gt_0.5_with_gwas_pval.tsv",
                        header=TRUE,sep='\t') %>%
  filter(gene_id == genei)
coloc_tiss = sort(unique(as.character(enloc_hits$tissue)))
coloc_gwas = sort(unique(as.character(enloc_hits$phenotype)))
write.table(enloc_hits, "APBB1IP_examp.coloc.txt",
            col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
lapply(unique(sort(as.character(enloc_hits$lead_snp))), function(vari) {
  LD(bed_matrix_gts,
     which(bed_matrix_gts@snps$id == overlap_var),
     which(bed_matrix_gts@snps$id == vari))
})

tabix_region=paste0(chri,":",posi-10^6,"-",posi+10^6)

## GWAS data
coloc_gwas
gwas_eos = tabix.read.table("input_files/tabixed/imputed_Astle_et_al_2016_Eosinophil_counts.txt.gz",
                      tabix_region) %>%
  mutate(gwas = "Astle_et_al_2016_Eosinophil_counts")
gwas_ret = tabix.read.table("input_files/tabixed/imputed_Astle_et_al_2016_High_light_scatter_reticulocyte_count.txt.gz",
                      tabix_region) %>%
  mutate(gwas = "Astle_et_al_2016_High_light_scatter_reticulocyte_count")
gwas_eos_baso = tabix.read.table("input_files/tabixed/imputed_Astle_et_al_2016_Sum_eosinophil_basophil_counts.txt.gz",
                           tabix_region) %>%
  mutate(gwas = "Astle_et_al_2016_Sum_eosinophil_basophil_counts")
gwas_lymp = tabix.read.table("input_files/tabixed/imputed_Astle_et_al_2016_Lymphocyte_counts.txt.gz",
                                 tabix_region) %>%
  mutate(gwas = "Astle_et_al_2016_Lymphocyte_counts")

all_gwas = Reduce('rbind', 
                  list(gwas_eos, gwas_ret, gwas_eos_baso, gwas_lymp)) %>%
  merge(all_ld,
        by.x=c('V2'),
        by.y=c('var1'))
all_gwas %>%
  arrange(r2) %>%
  ggplot(aes(V4,-log10(V11))) +
  geom_point(aes(fill=r2_cat,
                 shape=r2_cat=='lead'),
             size=1.2,
             stroke=.2) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(21,23)) +
  geom_label_repel(aes(label=ifelse(V1%in% c(overlap_var,sig_tiss_var), 
                                    as.character(var1), '')),
                   nudge_x=.1) +
  facet_wrap(~gwas,
              scales = 'free_y',
             nrow=4) +
  xlim(posi-10^5,posi+10^5)


## eQTL data
sig_tiss
coloc_tiss
eqtl_art = tabix.read.table("input_files/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/Artery_Tibial.allpairs.txt.gz",
                            tabix_region) %>%
  filter(gene_id == genei) %>%
  mutate(tiss='Artery_Tibial')
eqtl_pit = tabix.read.table("input_files/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/Pituitary.allpairs.txt.gz",
                            tabix_region) %>%
  filter(gene_id == genei) %>%
  mutate(tiss='Pituitary')
eqtl_thy = tabix.read.table("input_files/GTEx_Analysis_v8_eQTL_all_associations_indexed/results/Thyroid.allpairs.txt.gz",
                           tabix_region)  %>%
  filter(gene_id == genei) %>%
  mutate(tiss='Thyroid')

all_eqtl = Reduce('rbind', 
                  list(eqtl_art,eqtl_pit,eqtl_thy)) %>%
  mutate(pos=as.numeric(as.character(pos))) %>%
  merge(all_ld,
        by.x=c('variant_id'),
        by.y=c('var1'))
all_eqtl %>%
  arrange(r2) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r2_cat,
                 shape=r2_cat=='lead')) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(21,23)) +
  geom_label_repel(aes(label=ifelse(variant_id == overlap_var, 'overlap',
                                    ifelse(variant_id == sig_tiss_var[1], 'Art/Pit',
                                           ifelse(variant_id == sig_tiss_var[3], 'Thy', '')))), 
                   nudge_x=.1) +
  facet_wrap(~tiss,
             scales = 'free_y',
             nrow=3) +
  xlim(posi-10^5,posi+10^5)


rbind(all_gwas %>% 
        rename(variant_id=V2, chr=V3, pos=V4, pval_nominal=V11, r2=r2, r2_cat=r2_cat,
               data_type=gwas) %>%
        select(variant_id, chr, pos, pval_nominal, r2, r2_cat, data_type),
      all_eqtl %>%
        rename(data_type=tiss) %>%
        select(variant_id, chr, pos, pval_nominal, r2, r2_cat, data_type)) %>%
  mutate(label=ifelse(variant_id %in% c(overlap_var, sig_tiss_var), 'special', '')) %>%
  write.table("gwas_eqtl_pvals.txt",
              col.names=TRUE, row.names=FALSE,
              sep='\t', quote=FALSE)

### Final Plot
rbind(all_gwas %>% 
        filter(gwas=='Astle_et_al_2016_Eosinophil_counts') %>%
        rename(variant_id=V2, chr=V3, pos=V4, pval_nominal=V11, r2=r2, r2_cat=r2_cat) %>%
        select(variant_id, chr, pos, pval_nominal, r2, r2_cat) %>%
        mutate(data='Eosinophil counts GWAS'),
      all_eqtl %>%
        filter(tiss=='Thyroid') %>%
        select(variant_id, chr, pos, pval_nominal, r2, r2_cat) %>%
        mutate(data='Thyroid APBB1IP eQTL')) %>%
  mutate(label=ifelse(variant_id %in% c(overlap_var, sig_tiss_var), 'special', '')) %>%
  arrange(r2) %>%
  arrange(label) %>%
  ggplot(aes(pos,-log10(pval_nominal))) +
  geom_point(aes(fill=r2_cat,
                 shape=r2_cat=='lead',
                 size=label),
             stroke=.1) +
  geom_point(aes(x=ifelse(label=='special',pos,NA),
                 fill=r2_cat,
                 shape=r2_cat=='lead',
                 size=label),
             stroke=.5) +
  theme_classic() +
  scale_fill_manual(values=pal_locuszoom(palette = c("default"), alpha = 1)(6)) +
  scale_shape_manual(values=c(21,23)) +
  scale_size_manual(values=c(1,2)) +
  facet_wrap(~data,
             scales = 'free_y',
             nrow=3) +
  xlim(posi-10^5,posi+10^5)
ggsave("APBB1IP_examp.pdf",
       height=4,width=6)
  



## plot TF-eQTLs
gtex_gts = tabix.read.table("../TFi-eQTL/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz",
                            paste0(strsplit(overlap_var,"_")[[1]][1],":",
                                   strsplit(overlap_var,"_")[[1]][2],"-",
                                   strsplit(overlap_var,"_")[[1]][2])) %>%
  select(-c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='indiv')

art_tf_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[1],".tpm.gct"),
                            header=TRUE, sep='\t',
                            skip = 2) %>%
  filter(Description == tfi) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

art_gene_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[1],".tpm.gct"),
                            header=TRUE, sep='\t',
                            skip = 2) %>%
  filter(Name == genei) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

art_data = merge(art_tf_expr,art_gene_expr,
                  by=c('indiv','sample')) %>%
  merge(gtex_gts, by='indiv') %>%
  mutate(gt = factor(V1,
                     levels=c('0/0','0/1','1/1')))

art_corr = art_data %>%
  group_by(gt) %>%
  summarize(n_samp = n(),
            lm_b0 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[1],
            lm_b1 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[2])

art_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt)) +
  geom_abline(data=art_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','deepskyblue4','deepskyblue1')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle(paste(sig_tiss_short[1], '\n', 
                descriptioni, genei, '\n', 
                overlap_var))



pit_tf_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[2],".tpm.gct"),
                         header=TRUE, sep='\t',
                         skip = 2) %>%
  filter(Description == tfi) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

pit_gene_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[2],".tpm.gct"),
                           header=TRUE, sep='\t',
                           skip = 2) %>%
  filter(Name == genei) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

pit_data = merge(pit_tf_expr,pit_gene_expr,
                 by=c('indiv','sample')) %>%
  merge(gtex_gts, by='indiv') %>%
  mutate(gt = factor(V1,
                     levels=c('0/0','0/1','1/1')))
pit_data %>% 
  select(-indiv,-sample,-V1) %>%
  rename(gene_expr = V1.x, tf_expr = V1.y) %>%
  write.table("~/projects/MANUSCRIPT/data_tables/GxE_GWAS/APBB1IP_pituitary_data.txt",
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

pit_corr = pit_data %>%
  group_by(gt) %>%
  summarize(n_samp = n(),
            lm_b0 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[1],
            lm_b1 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[2])
pit_corr %>%
  write.table("~/projects/MANUSCRIPT/data_tables/GxE_GWAS/APBB1IP_pituitary_corr.txt",
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
  

pit_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt)) +
  geom_abline(data=pit_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle(paste(sig_tiss_short[2], '\n', 
                descriptioni, genei, '\n', 
                overlap_var))

pit_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=pit_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)'))
ggsave("APBB1IP_examp_tfieqtlpit.pdf",
       height=3,width=4)



thy_tf_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[3],".tpm.gct"),
                         header=TRUE, sep='\t',
                         skip = 2) %>%
  filter(Description == tfi) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

thy_gene_expr = read.table(paste0("../TFi-eQTL/phenotypes/tpm_by_tiss/",sig_tiss_short[3],".tpm.gct"),
                           header=TRUE, sep='\t',
                           skip = 2) %>%
  filter(Name == genei) %>%
  select(-c(Name,Description)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample') %>%
  separate(col=sample,sep='[.]',into=c('one','two'),remove=FALSE) %>%
  unite(col='indiv',c(one,two),sep='.')

thy_data = merge(thy_tf_expr,thy_gene_expr,
                 by=c('indiv','sample')) %>%
  merge(gtex_gts, by='indiv') %>%
  mutate(gt = factor(V1,
                     levels=c('0/0','0/1','1/1')))
thy_data %>% 
  select(-indiv,-sample,-V1) %>%
  rename(gene_expr = V1.x, tf_expr = V1.y) %>%
  write.table("~/projects/MANUSCRIPT/data_tables/GxE_GWAS/APBB1IP_thyroid_data.txt",
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

thy_corr = thy_data %>%
  group_by(gt) %>%
  summarize(n_samp = n(),
            lm_b0 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[1],
            lm_b1 = lm(log10(V1.y) ~ log10(V1.x))$coefficients[2])
thy_corr %>%
  write.table("~/projects/MANUSCRIPT/data_tables/GxE_GWAS/APBB1IP_thyroid_corr.txt",
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

thy_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt)) +
  geom_abline(data=thy_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)')) +
  ggtitle(paste(sig_tiss_short[3], '\n', 
                descriptioni, genei, '\n', 
                overlap_var))

thy_data %>%
  ggplot(aes(log10(V1.x),log10(V1.y))) +
  geom_point(aes(col=gt),
             size=.5) +
  geom_abline(data=thy_corr, aes(slope=lm_b1,intercept=lm_b0,col=gt)) +
  theme_classic() +
  scale_color_manual(values=c('black','olivedrab4','olivedrab2')) +
  xlab(paste0('log10(',tfi,' TPM)')) +
  ylab(paste0('log10(',descriptioni,' TPM)'))
ggsave("APBB1IP_examp_tfieqtlthy.pdf",
       height=3,width=4)

##
tf_info = read.table("~/projects/crosstiss_tf_corrs/input_files/TFs.info.curated.txt",
                     header=TRUE,sep='\t')

chips = read.table("input_files/chips.list",
           fill=TRUE)
info = read.table("input_files/tf_chipseq_20200110/info_files/final_chipseqs.list",
                  header=FALSE,sep='\t') %>%
  separate(V1, c(NA,'file'), "/")

chip_info = merge(chips, info,
      by.x='V1',
      by.y='file')

motif_files = list.files("input_files/motif_overlap_vars/",
                         "*.txt.gz$",
                         full.names = TRUE)
motif_info = do.call('rbind', lapply(motif_files, function(fii) {
  tfi=strsplit(strsplit(fii,"/")[[1]][4],"[.]")[[1]][1]
  tabix.read.table(fii,"chr10:26430000-26450000") %>%
    filter(V8=="Y") %>%
    mutate(tf=tfi)
}))

both_info = rbind(chip_info %>% mutate(tf=V3.y,start=V3.x,stop=V4.x,cell=V2.y,type='chip') %>% 
                    select(tf,start,stop,cell,type),
                  motif_info %>% mutate(start=V3,stop=V3,cell=NA,type='motif') %>%
                    select(tf,start,stop,cell,type)) %>%
  filter(tf %in% as.character(tf_info$TF)) %>%
  arrange(start)
offset=26438309
both_info_region = both_info %>%
  filter((start > offset-51 & start < offset+51) |
           (stop < offset+51 & stop > offset-51) |
           (start < offset-51 & stop > offset+51)) %>%
  mutate(start = ifelse(start<offset-50,offset-50,start),
         stop = ifelse(stop>offset+50,offset+50,stop)) %>%
  mutate(tf=factor(tf,levels=unique(as.character(tf))),
         tf_factor = as.numeric(tf))

ggplot(both_info_region) + 
  scale_x_continuous(position='top',limits=c(offset-51,offset+51)) +
  geom_rect(data=subset(both_info_region,type=='chip'),
            aes(xmin=start,xmax=stop,ymin=tf_factor,ymax=tf_factor+0.5))  +
  geom_vline(xintercept=offset,color='red') +
  geom_rect(data=subset(both_info_region,type=='motif'),
            aes(xmin=start-0.5,xmax=stop+0.5,ymin=tf_factor+0.1,ymax=tf_factor+0.4),
            fill='white',color='blue') +
  scale_y_continuous(name='TF', breaks=seq(1.3,length(levels(both_info_region$tf))+.3), 
                     labels=levels(both_info_region$tf)) +
  theme_classic() +
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank())


offset=26446058
both_info_region = both_info %>%
  filter((start > offset-51 & start < offset+51) |
           (stop < offset+51 & stop > offset-51) |
           (start < offset-51 & stop > offset+51)) %>%
  mutate(start = ifelse(start<offset-50,offset-50,start),
         stop = ifelse(stop>offset+50,offset+50,stop)) %>%
  mutate(tf=factor(tf,levels=unique(as.character(tf))),
         tf_factor = as.numeric(tf))

ggplot(both_info_region) + 
  scale_x_continuous(position='top',limits=c(offset-51,offset+51)) +
  geom_rect(data=subset(both_info_region,type=='chip'),
            aes(xmin=start,xmax=stop,ymin=tf_factor,ymax=tf_factor+0.5))  +
  geom_vline(xintercept=offset,color='red') +
  geom_rect(data=subset(both_info_region,type=='motif'),
            aes(xmin=start-0.5,xmax=stop+0.5,ymin=tf_factor+0.1,ymax=tf_factor+0.4),
            fill='white',color='blue') +
  scale_y_continuous(name='TF', breaks=seq(1.3,length(levels(both_info_region$tf))+.3), 
                     labels=levels(both_info_region$tf)) +
  theme_classic() +
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank())

