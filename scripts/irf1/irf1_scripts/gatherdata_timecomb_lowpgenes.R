#!/usr/bin/Rscript
##
##  gatherdata_timecomb_lowpgenes.R
##
##  EDF 5/11/21
##

library(dplyr)
library(ggplot2)
library(tidyr)
library(metap)

library(gaston)
library(seqminer)
library(tibble)
library(ggsci)
library(scales)
library(ggrepel)

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



moc_tf_eqtls = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                          header=TRUE,sep='\t')
eqtl_snps = read.table("eqtls/gene_var_gt.multi_or_cross.hets.txt",
                       header=FALSE, sep='\t') %>%
  mutate(chr_pos = paste(V3,V4,sep='_'))
IRF1_top_snps = read.table("input_files/TFi_hits_topvars.IRF1.bed",
                           header=FALSE,sep='\t') %>%
  mutate(chr_pos=paste(V1,V3,sep = '_'))
top_snps_IRF1_het_in_HEK = eqtl_snps %>% filter(chr_pos %in% IRF1_top_snps$chr_pos)
coding_snps = read.table("coding_snps/all_gene_exon_var_gt.hets.txt",
                         header=TRUE,sep = '\t')
gene_info = read.table("../crosstiss_tf_corrs/input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                       header=TRUE,sep='\t')


all_mpileups_indiv_annot = all_mpileups_indiv %>%
  mutate(totc=refc+altc,
         refp=refc/totc,
         altp=altc/totc) %>%
  merge(coding_snps %>% select(gene,exon,variant,chr,pos),
        by=c('chr','pos'), all.x=TRUE) %>%
  merge(gene_info %>% select(gene,description), by='gene') %>%
  mutate(time=factor(time,levels=c('0h','90m','12h')))



all_mpileups_comb_annot = all_mpileups_indiv_annot %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt,time,cond) %>%
  summarize(.groups='drop',
            num_samps = n(),
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





mpileups_comb_fi_test_lowp = mpileups_comb_fi_test %>%
  filter(fi_p < 0.05) %>%
  pivot_longer(cols=c(ref_C,alt_C,ref_P,alt_P)) %>%
  separate(name,sep='_',into=c('allele','cond'),remove=FALSE) %>%
  mutate(allele=factor(allele,levels=c('ref','alt')),
         cond=factor(cond,levels=c('P','C'))) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = max(value)*1.05) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  mutate(comb_p = ifelse(length(unique(fi_p)) > 1,
                         sumlog(unique(fi_p))$p,
                         first(fi_p)),
         gene_p_lab = paste(description,"p =",
                            ifelse(comb_p < 0.01, 
                                   round(comb_p,3),
                                   round(comb_p,2))))





gene_info = read.table("../crosstiss_tf_corrs/input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                       header=TRUE,sep='\t')

bed_matrix_gts = read.bed.matrix("/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/plink_genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01")

tiss_info = read.table("input_files/tissue_translation_colors_v8.txt",
                       header=TRUE,sep='\t')

irf1_overlap = read.table("../crosstiss_tf_corrs/overlaps/IRF1.txt.gz",
                          header=TRUE,sep='\t')



apply(unique(mpileups_comb_fi_test_lowp %>% ungroup() %>% select(gene, description) %>% 
               mutate(gene=as.character(gene), description=as.character(description))), 1, function(rowi) {
  
  genei = as.character(rowi[1])
  descriptioni = as.character(rowi[2])
  
  print(paste('Gene is',descriptioni,genei))
  
  cross_tiss_stats = read.table("../crosstiss_tf_corrs/correlations/by_tf/cross_tiss_tf_expr_corrs.med0.curated_set.MAF05.IRF1.top.adjp.txt",
                                  header=TRUE,sep='\t') %>%
    filter(gene == genei)
  
  sig_tiss_stats = do.call('rbind', lapply(list.files("../overlap/replication/indiv_tiss/",
                                                      "*fdr20*",
                                                      full.names = TRUE), 
                                           function(filei) {
                                             read.table(filei,
                                                        header=TRUE,sep='\t') %>%
                                               filter(tf=='IRF1',
                                                      phenotype_id == genei)
                                           } ))
  sig_tiss = sig_tiss_stats %>% pull(tiss) %>% as.character()
  sig_tiss_short = tiss_info %>%
    filter(TISSUE_NAME == sig_tiss[1]) %>%
    pull(TISSUE_ABBRV)
  
  cross_tiss_var = cross_tiss_stats %>% pull(var.min) %>% as.character()
  sig_tiss_var = sig_tiss_stats %>% pull(variant_id) %>% as.character()
  print(paste("Cross tiss top var is",cross_tiss_var,
              ', pvalue', cross_tiss_stats$sp_p.bf_meff_gene.min.bh_tf))
  print(paste(sig_tiss_stats %>% pull(tiss),"top var is",sig_tiss_var))
  
  if (length(cross_tiss_var) > 0 & length(sig_tiss_var) > 0) {
    print("Cross/within LD is")
    LD(bed_matrix_gts,
       which(bed_matrix_gts@snps$id == cross_tiss_var), 
       which(bed_matrix_gts@snps$id %in% sig_tiss_var)) %>%
      print()
  }
  
  
  irf1_overlap_gene = irf1_overlap %>% filter(gene==genei)
  
  irf1_overlap_both = irf1_overlap_gene %>% filter(tf_chip == 1, tf_motif ==1) %>% 
    pull(variant) %>% as.character()
  irf1_overlap_chip = irf1_overlap_gene %>% filter(tf_chip ==1) %>% 
    pull(variant) %>% as.character()
  irf1_overlap_motif = irf1_overlap_gene %>% filter(tf_motif ==1) %>% 
    pull(variant) %>% as.character()
  
  if (length(irf1_overlap_both) > 0) {
    print(paste('IRF1 overlap var(s) is/are',irf1_overlap_both))
    if (length(cross_tiss_var) > 0) {
      LD(bed_matrix_gts,
         which(bed_matrix_gts@snps$id == cross_tiss_var), 
         which(bed_matrix_gts@snps$id %in% irf1_overlap_both)) %>% print()
    }
    if (length(cross_tiss_var) > 0) {
      LD(bed_matrix_gts,
         which(bed_matrix_gts@snps$id %in% sig_tiss_var), 
         which(bed_matrix_gts@snps$id %in% irf1_overlap_both)) %>% print()
    }
  }
  
  print("any overlap LD")
  if (length(cross_tiss_var) > 0) {
    sapply(c(irf1_overlap_chip,irf1_overlap_motif), function(varj) {
      LD(bed_matrix_gts,
         which(bed_matrix_gts@snps$id == cross_tiss_var), 
         which(bed_matrix_gts@snps$id == varj)) %>% print()
    })
  }
  if (length(sig_tiss_var) > 0) {
    sapply(c(irf1_overlap_chip,irf1_overlap_motif), function(varj) {
      LD(bed_matrix_gts,
         which(bed_matrix_gts@snps$id %in% sig_tiss_var), 
         which(bed_matrix_gts@snps$id == varj)) %>% print()
    })
  }
  
})





