#!/usr/bin/Rscript
##
##  runFishersIRF1ASE_allgenes.R
##
##  EDF 4/15/21
##


library(dplyr)
library(ggplot2)
library(tidyr)
library(metap)

setwd("~/projects/IRF1/")

indiv_mpileups=list.files(path="coding_snps/", 
                          pattern="mpileup_.*_.*_.*[0-9]_allhetcodingsnps.txt.sum",
                          full.names = TRUE)
all_mpileups_indiv = do.call('rbind',lapply(indiv_mpileups, function(filei) {
  print(filei)
  read.table(filei,
             header=FALSE,sep='\t',quote='') %>%
    mutate(file=filei) %>%
    separate(file,sep='[/_]',c(NA,NA,NA,NA,'time','exp','samp')) %>%
    mutate(cond=gsub('[0-9]*','',exp))
}))

names(all_mpileups_indiv) <- c('chr','pos','ref','alt','refc','altc','otherc','reads','time','exp','samp','cond')


coding_snps = read.table("coding_snps/all_gene_exon_var_gt.hets.txt",
                         header=TRUE,sep = '\t')
length(unique(coding_snps$gene))


gene_info = read.table("../crosstiss_tf_corrs/input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                       header=TRUE,sep='\t')

irf1_expr = read.table("IRF1_norm_rc.txt",
                       header=FALSE,sep='\t') %>%
  separate(V1,c(NA,'samp'),sep='_',remove=FALSE) %>%
  mutate(log_irf1 = log2(V4))
names(irf1_expr) <- c('file','samp','tot_reads','irf1_rc','irf1_rc_norm','log_irf1')



all_mpileups_indiv_annot = all_mpileups_indiv %>%
  mutate(totc=refc+altc,
         refp=refc/totc,
         altp=altc/totc) %>%
  merge(coding_snps %>% select(gene,exon,variant,chr,pos),
        by=c('chr','pos'), all.x=TRUE) %>%
  merge(irf1_expr,by='samp') %>%
  merge(gene_info %>% select(gene,description), by='gene') %>%
  mutate(time=factor(time,levels=c('0h','90m','12h')))


all_mpileups_comb_annot = all_mpileups_indiv_annot %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt,time,cond) %>%
  summarize(.groups='drop',
            num_samps = n(),
            mean_irf1_rc_norm = mean(irf1_rc_norm),
            med_irf1_rc_norm = median(irf1_rc_norm),
            refc = sum(refc),
            altc = sum(altc),
            totc = sum(refc,altc),
            otherc = sum(otherc))


mpileups_comb_filt = all_mpileups_comb_annot %>%
  filter(cond %in% c('C','P')) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  filter('C' %in% cond, 'P' %in% cond) %>%
  group_by(gene,description,exon,variant) %>%
  filter(sum(totc) > 60,
         sum(otherc)/sum(totc) < 0.05,
         sum(altc)/sum(totc) > 0.05,
         sum(refc)/sum(totc) > 0.05)
  


mpileups_comb_fi_test = mpileups_comb_filt %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  summarize(.groups='drop',
            count=n(),
            ref_C = sum(refc[which(cond=='C')]),
            alt_C = sum(altc[which(cond=='C')]),
            tot_C = sum(totc[which(cond=='C')]),
            altp_C = alt_C/tot_C,
            ref_P = sum(refc[which(cond=='P')]),
            alt_P = sum(altc[which(cond=='P')]),
            tot_P = sum(totc[which(cond=='P')]),
            altp_P = alt_P/tot_P,
            fi_p = fisher.test(matrix(c(ref_C,alt_C,ref_P,alt_P),nrow=2))$p.value,
            chi_p = chisq.test(matrix(c(ref_C,alt_C,ref_P,alt_P),nrow=2))$p.value)

hist(mpileups_comb_fi_test$fi_p)
hist(mpileups_comb_fi_test$chi_p)

mpileups_comb_fi_test %>%
  filter(fi_p < 0.05)
mpileups_comb_fi_test %>%
  filter(chi_p < 0.05)

# 
# 
# mpileups_comb_fi_test_allps = mpileups_comb_fi_test %>%
#   group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
#   summarize(.groups='drop',
#             count=n(),
#             p_0h = if ('0h' %in% time) {fi_p[which(time=='0h')]} else {NA},
#             p_90m = if ('90m' %in% time) {fi_p[which(time=='90m')]} else {NA},
#             p_12h = if ('12h' %in% time) {fi_p[which(time=='12h')]} else {NA},
#             comb_p = ifelse(count>1,sumlog(c(fi_p))$p,fi_p),
#             min_p = min(fi_p)) %>%
#   mutate(comb_p_bh = p.adjust(comb_p, method='BH'))
# 
# hist(mpileups_comb_fi_test_allps$comb_p)
# hist(mpileups_comb_fi_test_allps$min_p)
# 
# hist(mpileups_comb_fi_test_allps$p_0h)
# hist(mpileups_comb_fi_test_allps$p_90m)
# hist(mpileups_comb_fi_test_allps$p_12h)



write.table(mpileups_comb_fi_test,
            "results/Fi_tst_combined_byvar_timecomb.allgenes.txt",
            col.names=TRUE, row.names=FALSE, sep='\t', quote = FALSE)


# mpileups_comb_fi_test_allps_bygene = mpileups_comb_fi_test %>%
#   group_by(gene,description) %>%
#   summarize(.groups='drop',
#             count=n(),
#             n_vars = length(unique(variant)),
#             min_p_0h = min(fi_p[which(time=='0h')]),
#             min_p_90m = min(fi_p[which(time=='90m')]),
#             min_p_12h = min(fi_p[which(time=='12h')]),
#             comb_p = ifelse(count>1,sumlog(c(fi_p))$p,fi_p),
#             min_p = min(fi_p),
#             min_p_exon = exon[which(fi_p==min_p)],
#             min_p_var = variant[which(fi_p==min_p)]) %>%
#   mutate(comb_p_bh = p.adjust(comb_p, method='BH'),
# #          comb_gt1_p_bh = p.adjust(ifelse(count>1,comb_p,NA),method='BH'))
# 
# 
# write.table(mpileups_comb_fi_test_allps_bygene,
#             "results/Fi_tst_combined_bygene.allgenes.txt",
#             col.names=TRUE, row.names=FALSE, sep='\t', quote = FALSE)
# 
# hist(mpileups_comb_fi_test_allps_bygene$comb_p)
# hist(mpileups_comb_fi_test_allps_bygene %>% filter(count>1) %>% pull(comb_p))
# 
# mpileups_comb_fi_test_allps_bygene %>%
#   filter(comb_p < 0.05)



# 
# 
# 
# moc_tf_eqtls = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
#                             header=TRUE,sep='\t')
# 
# table(unique(mpileups_comb_fi_test_allps$gene) %in% unique(mpileups_comb_fi_test_allps %>% filter(comb_p < 0.05) %>% pull(gene)),
#       unique(mpileups_comb_fi_test_allps$gene) %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
# 
# eqtl_snps = read.table("eqtls/gene_var_gt.multi_or_cross.hets.txt",
#                        header=FALSE, sep='\t') %>%
#   mutate(chr_pos = paste(V3,V4,sep='_'))
# top_snps = read.table("input_files/TFi_hits_topvars.IRF1.bed",
#                       header=FALSE,sep='\t') %>%
#   mutate(chr_pos=paste(V1,V3,sep = '_'))
# 
# top_snps_het_in_HEK = eqtl_snps %>% filter(chr_pos %in% top_snps$chr_pos)
# 
# table(unique(mpileups_comb_fi_test_allps$gene) %in% unique(mpileups_comb_fi_test_allps %>% filter(comb_p < 0.05) %>% pull(gene)),
#       unique(mpileups_comb_fi_test_allps$gene) %in% unique(top_snps_het_in_HEK$V1))
# 
# 
# mpileups_comb_fi_test_allps %>% filter(gene %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
# mpileups_comb_fi_test_allps %>% filter(gene %in% unique(top_snps_het_in_HEK$V1))
# 
# mpileups_comb_fi_test_allps_bygene %>% filter(gene %in% unique(moc_tf_eqtls %>% filter(tf == 'IRF1') %>% pull(phenotype_id)))
# mpileups_comb_fi_test_allps_bygene %>% filter(gene %in% unique(top_snps_het_in_HEK$V1))

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



# irf1_reg = read.table("input_files/merged_interactions.txt",
#                       header=TRUE,sep='\t') %>%
#   filter(TF=='IRF1')
# 
# irf_reg = read.table("input_files/merged_interactions.txt",
#                      header=TRUE,sep='\t') %>%
#   filter(startsWith(as.character(TF),'IRF'))
# 
# stat_reg = read.table("input_files/merged_interactions.txt",
#                       header=TRUE,sep='\t') %>%
#   filter(startsWith(as.character(TF),'STAT'))
# 
# nfkb_reg = read.table("input_files/merged_interactions.txt",
#                       header=TRUE,sep='\t') %>%
#   filter(startsWith(as.character(TF),'NFKB'))
