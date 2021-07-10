#!/usr/bin/Rscript
##
##  plotIRF1ASE_timecomb_lowpgenes.R
##
##  EDF 4/20/21
##


library(dplyr)
library(ggplot2)
library(tidyr)
library(metap)
library(UpSetR)

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
IRF1_top_snps = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.vars_unique.txt",
                           header=TRUE,sep='\t') %>%
  filter(tf=='IRF1') %>%
  separate(variant,into=c('chr','pos'),sep='_',remove=FALSE) %>%
  unite(chr_pos,chr,pos)
# 
# IRF1_top_snps = read.table("input_files/TFi_hits_topvars.IRF1.bed",
#                            header=FALSE,sep='\t') %>%
#   mutate(chr_pos=paste(V1,V3,sep = '_'))
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
hist(mpileups_comb_fi_test$fi_p,
     main='Histogram of IRF1 knockdown ASE vs. condition p values',
     xlab="Fisher's exact test p value",
     ylab='Genes')



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
write.table(mpileups_comb_fi_test_lowp,
            "results/Fi_tst_combined_byvar_timecomb.lowpgenes.txt",
            col.names = TRUE, row.names = FALSE, sep='\t', quote=FALSE)


mpileups_comb_fi_test_lowp %>%
  ggplot(aes(1,value)) +
  geom_col(aes(alpha=cond, fill=allele),
           position=position_dodge(),
           width=.75) +
  geom_text(aes(1,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  scale_fill_manual(values=c('coral','cornflowerblue')) +
  geom_text(aes(1,p_height*1.1,
                label=NA)) +
  facet_wrap(~gene_p_lab,
             scales = 'free_y') +
  theme_classic() +
  ylab('Read Count') +
  xlab('Time Point')








mpileups_comb_filt_times = all_mpileups_comb_annot %>%
  filter(cond %in% c('C','P')) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt,time) %>%
  filter('C' %in% cond, 'P' %in% cond) %>%
  group_by(gene,description,exon,variant) %>%
  filter(sum(totc) > 60,
         sum(otherc)/sum(totc) < 0.05,
         sum(altc)/sum(totc) > 0.05,
         sum(refc)/sum(totc) > 0.05)


mpileups_comb_fi_test_times = mpileups_comb_filt_times %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt,time) %>%
  summarize(.groups='drop',
            count=n(),
            ref_C = refc[which(cond=='C')],
            alt_C = altc[which(cond=='C')],
            tot_C = totc[which(cond=='C')],
            altp_C = alt_C/tot_C,
            ref_P = refc[which(cond=='P')],
            alt_P = altc[which(cond=='P')],
            tot_P = totc[which(cond=='P')],
            altp_P = alt_P/tot_P,
            fi_p = fisher.test(matrix(c(ref_C,alt_C,ref_P,alt_P),nrow=2))$p.value)





mpileups_comb_fi_test_times2 = mpileups_comb_fi_test_times %>%
  pivot_longer(cols=c(ref_C,alt_C,ref_P,alt_P)) %>%
  separate(name,sep='_',into=c('allele','cond'),remove=FALSE) %>%
  mutate(allele=factor(allele,levels=c('ref','alt')),
         cond=factor(cond,levels=c('P','C'))) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt,time) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = max(value)*1.05) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  mutate(comb_p = ifelse(length(unique(fi_p)) > 1,
                         sumlog(unique(fi_p))$p,
                         first(fi_p)),
         gene_p_lab = paste(description,"p =",round(comb_p,2)))



mpileups_comb_fi_test_lowp_plot = merge(mpileups_comb_fi_test_times2 %>% select(-c(comb_p,gene_p_lab)),
                                        mpileups_comb_fi_test_lowp %>% select(gene,description,exon,variant,allele,chr,pos,ref,alt,
                                                                              cond,fi_p,comb_p,gene_p_lab),
                                        by=c('gene','description','exon','variant','chr','pos','ref','alt','allele','cond'))


mpileups_comb_fi_test_lowp_plot %>%
  ggplot(aes(time,value)) +
  geom_col(aes(alpha=cond, fill=allele),
           position=position_dodge(),
           width=.75) +
  geom_text(aes(time,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  scale_fill_manual(values=c('coral','cornflowerblue')) +
  geom_text(aes(time,p_height*1.1,
                label=NA)) +
  facet_wrap(~gene_p_lab,
             scales = 'free_y') +
  theme_classic() +
  ylab('Read Count') +
  xlab('Time Point')

mpileups_comb_fi_test_lowp_plot %>%
  pivot_longer(cols=c('tot_C','altp_C','tot_P','altp_P'),names_to='cat',values_to='val') %>%
  separate(cat,into=c('var','cond2')) %>%
  filter(cond==cond2) %>%
  pivot_wider(names_from='var',values_from='val') %>%
  ggplot(aes(time,altp)) +
  geom_point(aes(color=cond,group=cond,size=tot),
             position=position_dodge(width=.5)) +
  geom_line(aes(color=cond,group=cond),
            position=position_dodge(width=.5)) +
  geom_text(aes(time,1,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  scale_color_manual(values=c('coral','cornflowerblue')) +
  geom_text(aes(time,1,
                label=NA)) +
  facet_wrap(~gene_p_lab) +
  ylim(0,1) +
  theme_classic() +
  ylab('Alt Allelic Fraction') +
  xlab('Time Point')

mpileups_comb_fi_test_lowp_plot %>%
  pivot_longer(cols=c('tot_C','altp_C','tot_P','altp_P'),names_to='cat',values_to='val') %>%
  separate(cat,into=c('var','cond2')) %>%
  filter(cond==cond2) %>%
  pivot_wider(names_from='var',values_from='val') %>%
  mutate(afc = (altp) / (1-altp),
         log2afc = log2(afc)) %>%
  group_by(time,gene_p_lab) %>%
  mutate(pheightafc = mean(log2afc)) %>%
  ggplot(aes(time,log2afc)) +
  geom_point(aes(color=cond,group=cond,size=tot)) +
  geom_line(aes(color=cond,group=cond)) +
  geom_text(aes(time,pheightafc,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  scale_color_manual(values=c('coral','cornflowerblue')) +
  geom_text(aes(time,0,
                label=NA)) +
  facet_wrap(~gene_p_lab,
             scales = 'free_y') +
  theme_classic() +
  ylab('Log2(Allelic Fold Change)') +
  xlab('Time Point')




tested_genes = as.character(unique(mpileups_comb_fi_test$gene))
HEK_lowp_genes = as.character(unique(mpileups_comb_fi_test_lowp$gene))

IRF1_hits_cross = read.table("../overlap/replication/sig_assoc.fdr05.cross_hits.txt",
                             header=TRUE,sep='\t') %>%
  filter(phenotype_id %in% tested_genes,
         tf=='IRF1')

IRF1_hits_within = read.table("../overlap/replication/sig_assoc.fdr20.any_tiss_hits.txt",
                              header=TRUE, sep='\t' ) %>%
  filter(phenotype_id %in% tested_genes,
         tf=='IRF1')

pdf("scripts/upset_20210604.pdf",width=4,height=4)
upset(fromList(list(
  cross = IRF1_hits_cross %>% 
    pull(phenotype_id) %>% as.character(),
  tiss = IRF1_hits_within %>%
    pull(phenotype_id) %>% as.character(),
  HEK = HEK_lowp_genes)))
dev.off()





# 
# temp_asb = read.table("../ASB/ADASTRA_data/ADASTRA_Susan/release_dump/TF/IRF1_HUMAN.tsv",
#                       header=FALSE,sep='\t') %>%
#   unite(chr_pos,V1,V2,remove=FALSE)
# temp_asb[which(temp_asb$chr_pos %in% IRF1_top_snps$chr_pos),]
# temp_asb[which(temp_asb$chr_pos %in% eqtl_snps$chr_pos),]
# temp_asb[which(temp_asb$V15 < 0.1 | temp_asb$V16 < 0.1),]
# 
# temp_asb2 = read.table("../ASB/ADASTRA_Susan/by_tf/IRF1_overview.txt",
#                        header=TRUE, sep='\t') %>%
#   separate(variant,c('chr','pos'),'_') %>%
#   unite(chr_pos,chr,pos) %>%
#   filter(has_ASB_fdr10 > 0)



