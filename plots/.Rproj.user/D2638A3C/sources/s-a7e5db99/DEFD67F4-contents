#!/usr/bin/Rscript
##
##  plotASB_Susan_20210511.R
##
##  EDF 5/11/2021
##

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

setwd("~/projects/ASB/")


tf_info = read.table("input_files/TFs.info.curated.txt",
                     header=TRUE,
                     sep='\t')

cor_genes = read.table("../overlap/replication/sig_assoc.fdr20.multi_or_cross_hits.txt",
                       header=TRUE,sep='\t') %>%
  mutate(tf_gene = paste(tf,phenotype_id))
gene_var = read.table("input_files/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                      header=TRUE,sep='\t') %>%
  separate(var,c(NA,NA,'ref','alt'), remove=FALSE)

#you can get enrichment per TF from 
#/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/ADASTRA/ADASTRA_data.txt. 
#if the TF-eQTL pair has ASB and the TFs match, then tf_match = 1 and has_ASB = 1; the TF is under the accB_tf column

#the combined permutation data for the TF-matched analysis is in 
#/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/perm_stats/perms_combined_match.txt. 
#i’ve shown perm plots for diff and avg_corr (which is just diff without the uncorr part aka what we talked about with harmen)

#the permutation data for the not-matched analysis is in 
#/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/combined_ASB/perms_combined_overall.txt. 
#i plotted eq8_ASB for the corr genes (sig_genes = TRUE)



### ASB per TF
asb_tf = read.table("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/ADASTRA_Susan_new/perm_stats_match/ADASTRA_Susan_combined_eq2.txt",
                    header=TRUE,sep='\t')

p_perm_match = (sum(asb_tf$avg_corr_05>(asb_tf%>%filter(perm==0)%>%pull(avg_corr_05)))+1)/5000

asb_tf %>%
  mutate(temp='temp') %>%
  pivot_wider(id_cols=temp, names_from = perm, values_from=avg_corr_05) %>%
  pivot_longer(cols=as.character(seq(1,10000)),
               names_to='perm', values_to='avg_corr') %>%
  rename(real=`0`) %>%
  ggplot(aes(avg_corr)) +
  geom_histogram(fill='gray60') +
  geom_vline(aes(xintercept=real)) +
  theme_classic() +
  xlab('ASB Enrichment') +
  annotate(geom='text',x=-0.015,y=850,
           label=paste('perm p =',p_perm_match),
           hjust=0) +
  annotate(geom='text',
           x=asb_tf%>%filter(perm==0)%>%pull(avg_corr_05)+0.0005,
           y=950,
           label=round(asb_tf%>%filter(perm==0)%>%pull(avg_corr_05),3),
           hjust=0)




asb_all = read.table("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/ADASTRA_Susan_new/combined_ASB/combined.perm_stats.txt",
                     header=TRUE,sep='\t') %>%
  filter(sig_genes == TRUE)


p_perm_all = (sum(asb_all$eq2_ASB_05>(asb_all%>%filter(perm==0)%>%pull(eq2_ASB_05)))+1)/5000

asb_all %>%
  mutate(temp='temp') %>%
  pivot_wider(id_cols=temp, names_from = perm, values_from=eq2_ASB_05) %>%
  pivot_longer(cols=as.character(seq(1,10000)),
               names_to='perm', values_to='avg_corr') %>%
  rename(real=`0`) %>%
  ggplot(aes(avg_corr)) +
  geom_histogram(fill='gray60') +
  geom_vline(aes(xintercept=real)) +
  theme_classic() +
  xlab('ASB Enrichment') +
  annotate(geom='text',x=-0.06,y=850,
           label=paste('perm p =',p_perm_all),
           hjust=0) +
  annotate(geom='text',
           x=asb_all%>%filter(perm==0)%>%pull(eq2_ASB_05)+0.002,
           y=950,
           label=round(asb_all%>%filter(perm==0)%>%pull(eq2_ASB_05),3),
           hjust=0)




# asb_bytf = read.table("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/ADASTRA_Susan/perm_stats_match/ADASTRA_Susan_combined_eq2.txt",
#                       header=TRUE,sep='\t') %>%
#   mutate(run='bytf') %>%
#   pivot_wider(id_cols=run, names_from = perm, values_from=avg_corr_05) %>%
#   pivot_longer(cols=as.character(seq(1,10000)),
#                names_to='perm', values_to='eq2_avg') %>%
#   rename(real=`0`) %>%
#   group_by(run) %>%
#   summarize(real=first(real),
#             num_gt = sum(eq2_avg >= real),
#             p_perm = (num_gt+1)/10001*2)
# 
# asb_together = read.table("/gpfs/commons/groups/lappalainen_lab/atsu/projects/ASB/ADASTRA/combined_ASB/perms_combined_overall.txt",
#                           header=TRUE,sep='\t') %>%
#   pivot_wider(id_cols=perm, names_from=sig_genes, values_from=eq8_ASB) %>%
#   mutate(diff=`TRUE`-`FALSE`,
#          run='together') %>%
#   pivot_wider(id_cols=run, names_from=perm, values_from=diff) %>%
#   pivot_longer(cols=as.character(seq(1,10000)),
#                names_to='perm', values_to='diff') %>%
#   rename(real=`0`) %>%
#   group_by(run) %>%
#   summarize(real=first(real),
#             num_gt = sum(diff >= real),
#             p_perm = (num_gt+1)/10001*2)
# 
# 
# rbind(asb_together,
#       asb_bytf) %>%
#   mutate(p_lab = paste('p =',round(p_perm,3))) %>%
#   mutate(run = factor(run, levels=c('together','bytf'))) %>%
#   ggplot(aes(run,real)) +
#   geom_col(width=.02) +
#   geom_point(col='darkorchid4',
#              size=2) +
#   geom_hline(yintercept=0) +
#   geom_text(aes(run,
#                 real+0.005,
#                 label=p_lab),
#             size=3) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=15,hjust=1),
#         strip.text = element_text(hjust=0)) +
#   ylab("Overlap Enrichment")
# ggsave('../overlap/figs3/ASB_overlap.20210413.pdf',
#        height=2,width=2)


tf_files = list.files("ADASTRA_Susan/by_tf",
                      full.names = TRUE)

asb_by_tf = do.call('rbind', lapply(tf_files, function(filei) {
  
  tfi = strsplit(strsplit(filei,'/')[[1]][3],'[_]')[[1]][1]
  
  tablei = read.table(filei, header=TRUE, sep='\t')
  
  p_ASB = sum(tablei$has_ASB_fdr05)/sum(!is.na(tablei$has_ASB_fdr05))
  
  tablei %>%
    group_by(gene,corr_gene,corr_tf) %>%
    summarise(tf = tfi,
              exp = p_ASB*n(),
              obs = sum(has_ASB_fdr05),
              diff = obs-exp)
}) )

asb_by_tf %>%
  filter(corr_gene == 1,
         grepl(tf, corr_tf)) %>%
  group_by(tf) %>%
  summarize(sum_diff = sum(diff),
            num_corr = sum(!is.na(diff)),
            avg_diff = sum_diff/num_corr)




tf_perm_files = list.files("ADASTRA_Susan_new/perm_stats_match",
                           pattern="*.perm_stats.txt",
                      full.names = TRUE)

asb_perm_by_tf = do.call('rbind', lapply(tf_perm_files, function(filei) {
  
  tfi = strsplit(strsplit(filei,'/')[[1]][3],'[.]')[[1]][1]
  print(tfi)
  
  tablei = read.table(filei, header=TRUE, sep='\t')
  
  if(sum(tablei$sig_genes) > 0) {
    tablei %>%
    filter(sig_genes) %>%
    pivot_wider(id_cols = num_sig, 
                names_from = perm,
                values_from = eq2_ASB_05) %>%
    pivot_longer(cols=as.character(seq(1,10000)),
                 names_to='perm',
                 values_to='perm_stat') %>%
    rename(real_stat='0') %>%
    mutate(tf = tfi)
  } 
}) )

asb_perm_by_tf_sum = asb_perm_by_tf %>%
  group_by(tf, real_stat) %>%
  summarize(num_gte = sum(perm_stat >= real_stat),
    p_perm_tf = ( min( ( min( sum(perm_stat >= real_stat),
                               sum(perm_stat <= real_stat) ) + 1 )/ 5000, 1) ))

asb_perm_by_tf %>%
  ggplot() +
  geom_histogram(aes(perm_stat)) +
  geom_vline(aes(xintercept=real_stat)) +
  facet_wrap(~tf,
             scales = 'free') +
  theme_classic()

asb_perm_by_tf_sum %>%
  ggplot(aes(fct_reorder(tf, real_stat),real_stat)) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  geom_col(width=.2) +
  scale_color_gradient(low='black',high='darkorchid2') +
  coord_flip() +
  theme_classic()




ADASTRA_susan_table = do.call('rbind', apply(tf_info, 1, function(rowi) {
  tfi=rowi[1]
  print(tfi)
  filei=list.files("ADASTRA_data/ADASTRA_Susan/release_dump/TF/", paste0("^",tfi,"_HUMAN.tsv"), full.names = TRUE)
  if (length(filei) > 0) {
    read.table(filei,
               header=TRUE,comment.char = '',sep='\t') %>%
      mutate(tf=tfi)
  }
}) )


ADASTRA_susan_annot = merge(gene_var,
                            ADASTRA_susan_table %>% select(X.chr,pos,ID,ref,alt,fdrp_bh_ref,fdrp_bh_alt,tf),
                            by.x=c('chr','pos','ref','alt'),
                            by.y=c('X.chr','pos','ref','alt')) %>% 
  mutate(ASB_either = fdrp_bh_ref < 0.05 | fdrp_bh_alt < 0.05,
         corr = paste(tf,gene) %in% cor_genes$tf_gene)


fisher.test(table(ADASTRA_susan_annot$ASB_either,ADASTRA_susan_annot$corr))



ADASTRA_susan_bygene = ADASTRA_susan_annot %>%
  group_by(tf,gene) %>%
  summarize(n_var = n(),
            any_ASB = any(ASB_either),
            corr = first(corr))

fisher.test(table(ADASTRA_susan_bygene$any_ASB,ADASTRA_susan_bygene$corr))


ADASTRA_susan_tfsum = ADASTRA_susan_bygene %>%
  group_by(tf) %>%
  summarize(n_gene = n(),
            n_ASB = sum(any_ASB),
            n_corr = sum(corr),
            exp_ASB_corr = n_ASB*n_corr/n_gene,
            n_ASB_corr = sum(any_ASB & corr),
            fi_OR = fisher.test(table(factor(any_ASB, levels=c(FALSE,TRUE)),
                                      factor(corr, levels=c(FALSE,TRUE))))$estimate,
            fi_logOR = log2(fi_OR),
            fi_p = fisher.test(table(factor(any_ASB, levels=c(FALSE,TRUE)),
                                     factor(corr, levels=c(FALSE,TRUE))))$p.value)

ADASTRA_susan_tfsum %>%
  mutate(tf = fct_reorder(tf, fi_logOR)) %>%
  filter(!is.na(fi_logOR)) %>%
  ggplot(aes(tf,fi_logOR)) +
  geom_col(width=.1) +
  geom_point(aes(col=fi_p < 0.05,size=exp_ASB_corr)) +
  theme_classic() +
  scale_color_manual(values=c('gray','darkorchid4')) +
  coord_flip()



asb_annot = merge(asb_perm_by_tf_sum, ADASTRA_susan_tfsum,
      by='tf')

asb_annot %>%
  filter(exp_ASB_corr >= 2) %>%
  ggplot(aes(fct_reorder(tf, real_stat),real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient(low='gray',high='darkorchid4') +
  coord_flip() +
  theme_classic()

asb_annot %>%
  filter(exp_ASB_corr < 2) %>%
  summarize(new_real_stat = sum(real_stat * n_corr) / sum(n_corr))



asb_perm_lt2 = asb_perm_by_tf %>%
  filter(tf %in% as.character(asb_annot %>% filter(exp_ASB_corr < 2) %>% pull(tf))) %>%
  merge(ADASTRA_susan_tfsum,
        by='tf') %>%
  group_by(perm) %>%
  summarize(new_perm_stat = sum(perm_stat * num_sig) / sum(num_sig),
            new_real_stat = sum(real_stat * num_sig / sum(num_sig)),
            new_exp_ASB_corr = sum(exp_ASB_corr)) %>%
  ungroup() %>%
  summarize(tf = 'other',
            real_stat = unique(new_real_stat),
            exp_ASB_corr = unique(new_exp_ASB_corr),
            num_gte = sum(new_perm_stat >= new_real_stat),
            p_perm_tf = ( min( ( min( sum(new_perm_stat >= new_real_stat),
                                      sum(new_perm_stat <= new_real_stat) ) + 1 )/ 5000, 1) ))


asb_perm_lt1 = asb_perm_by_tf %>%
  filter(tf %in% as.character(asb_annot %>% filter(exp_ASB_corr < 1) %>% pull(tf))) %>%
  merge(ADASTRA_susan_tfsum,
        by='tf') %>%
  group_by(perm) %>%
  summarize(new_perm_stat = sum(perm_stat * num_sig) / sum(num_sig),
            new_real_stat = sum(real_stat * num_sig / sum(num_sig)),
            new_exp_ASB_corr = sum(exp_ASB_corr)) %>%
  ungroup() %>%
  summarize(tf = 'other',
            real_stat = unique(new_real_stat),
            exp_ASB_corr = unique(new_exp_ASB_corr),
            num_gte = sum(new_perm_stat >= new_real_stat),
            p_perm_tf = ( min( ( min( sum(new_perm_stat >= new_real_stat),
                                      sum(new_perm_stat <= new_real_stat) ) + 1 )/ 5000, 1) ))


rbind(asb_annot %>%
        filter(exp_ASB_corr >= 1) %>%
        select(tf,real_stat,exp_ASB_corr,num_gte,p_perm_tf),
      asb_perm_lt1) %>%
  ggplot(aes(fct_reorder(tf, real_stat),real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient(low='gray',high='darkorchid4',
                       limits=c(0,2)) +
  coord_flip() +
  theme_classic()

rbind(asb_annot %>%
        filter(exp_ASB_corr >= 1) %>%
        select(tf,real_stat,exp_ASB_corr,num_gte,p_perm_tf),
      asb_perm_lt1) %>%
  write.table("ASB_susan_filtexpASB1.txt",
              col.names=TRUE,row.names=FALSE,
              sep='\t',quote=FALSE)

rbind(asb_annot %>%
        filter(exp_ASB_corr >= 2) %>%
        select(tf,real_stat,exp_ASB_corr,num_gte,p_perm_tf),
      asb_perm_lt2) %>%
  ggplot(aes(fct_reorder(tf, real_stat),real_stat,
             size=exp_ASB_corr)) +
  geom_col(width=.05) +
  geom_point(aes(col = -log10(p_perm_tf))) +
  scale_color_gradient(low='gray',high='darkorchid4',
                       limits=c(0,2)) +
  coord_flip() +
  theme_classic()

rbind(asb_annot %>%
        filter(exp_ASB_corr >= 2) %>%
        select(tf,real_stat,exp_ASB_corr,num_gte,p_perm_tf),
      asb_perm_lt2) %>%
  write.table("ASB_susan_filtexpASB2.txt",
              col.names=TRUE,row.names=FALSE,
              sep='\t',quote=FALSE)




ADASTRA_susan_byvar_any = ADASTRA_susan_annot %>%
  group_by(var) %>%
  summarize(n_var = n(),
            any_ASB = any(ASB_either))
table(ADASTRA_susan_byvar_any$any_ASB)

data.frame(tf='all',
           real_stat=asb_all%>%filter(perm==0)%>%pull(eq2_ASB_05),
           exp_ASB_corr = asb_all[1,]$num_sig * 
             sum(ADASTRA_susan_byvar_any$any_ASB) / nrow(ADASTRA_susan_byvar_any),
           num_gte = sum(asb_all$eq2_ASB_05>(asb_all%>%filter(perm==0)%>%pull(eq2_ASB_05))),
           p_perm_tf = p_perm_all ) %>%
  write.table("ASB_susan_allcombined.txt",
              col.names=TRUE,row.names=FALSE,
              sep='\t',quote=FALSE)
  
  
  
