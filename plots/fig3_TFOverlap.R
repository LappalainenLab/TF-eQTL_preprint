#!/usr/bin/Rscript
##
##  fig3_TFOverlap.R
##
##  EDF 5/7/21
##


library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)


setwd("~/projects/MANUSCRIPT/")

perm_files=list.files("data_tables/overlap","perms.*",
                      full.names = TRUE)


overlap_sum = do.call('rbind', lapply(perm_files, function(filei) {
  datai = read.table(filei,
                     header=TRUE,sep='\t') %>%
    pivot_longer(cols=c(true_stat,diff_stat,diff_stat_2),
                 names_to="var", values_to="value") %>%
    pivot_wider(id_cols=c(name,var,count_true,count_diff),
                names_from=perm,
                values_from=value) %>%
    pivot_longer(cols=as.character(1:10^4), 
                 names_to = c('perm'), 
                 values_to = 'stat') %>%
    rename(real = `0`)
  
  datai %>%
    group_by(name,var,real) %>%
    summarize(.groups='drop',
              tf_count = unique(ifelse(var=="true_stat",count_true,
                                       ifelse(var=="diff_stat",count_diff,
                                              ifelse(var=='diff2_stat',count_diff,NA)))),
              perms_lt = sum(stat <= real),
              perms_gt = sum(stat >= real),
              perm_p = ( min(perms_lt,perms_gt) + 1 ) / 5000.5,
              perm_mean = mean(stat),
              perm_sd = sqrt(var(stat)),
              perm_z = ( unique(real) - perm_mean ) / perm_sd ) %>%
    separate(name,c("eq",'data_type')) %>%
    mutate(file = strsplit(filei,'[/]')[[1]][3],
           data = paste(eq,data_type,var),
           data_type = factor(data_type,
                              levels = c('chip','motif','both')))
}))

overlap_sum %>%
  separate(file,into=c(NA,NA,NA,"list",NA),sep='[.]',remove = FALSE) %>%
  mutate(list = factor(list, 
                       levels=c('cross_hits','protein_hits_filtered','any_tiss_hits','tiss_cross_hits','multi_tiss_hits','multi_or_cross_hits')),
         list_lab = factor(ifelse(list == 'cross_hits', 'Cross\nExpr',
                                  ifelse(list=='protein_hits_filtered','Cross\nProt',
                                         ifelse(list=='any_tiss_hits','\nWithin\n',
                                                ifelse(list=='multi_or_cross_hits','2+ Within\nor\nWithin + Cross','')))),
                           levels=c('Cross\nExpr','Cross\nProt','\nWithin\n','2+ Within\nor\nWithin + Cross')),
         data_type_lab = factor(ifelse(data_type == 'chip', 'ChIPseq Overlap',
                                       ifelse(data_type == 'motif','Motif Overlap',
                                              ifelse(data_type == 'both','Both Overlap',''))),
                                levels=c('ChIPseq Overlap','Motif Overlap','Both Overlap'))) %>%
  filter(!is.na(list)) %>%
  mutate(p_lab = ifelse(perm_p < 0.0002, "p < 2e-04",
                        ifelse(perm_p < 0.001, paste('p =',format(round(perm_p,4),scientific = FALSE)),
                               ifelse(perm_p < 0.005, paste('p =',format(round(perm_p,3),scientific = FALSE)),
                                      paste('p =',round(perm_p,2))))) ) %>%
  filter(eq=='eq2',
         var == 'true_stat') %>%
  ggplot(aes(list, real)) +
  geom_col(width=.02) +
  geom_hline(yintercept=0) +
  geom_point(col='darkorchid4') +
  geom_point(aes(list,ifelse(real > 0, real*1.5, real)),col=NA) +
  geom_text(aes(list,
                rep(c(0.002,0.09,0.006),6),
                label=p_lab),
            size=3) +
  facet_wrap(~data_type_lab, scales='free_y', nrow=3) +
  theme_classic() +
  theme(strip.text = element_text(hjust=0)) +
  ylab("Overlap Enrichment") +
  xlab('Dataset')


overlap_sum %>%
  separate(file,into=c(NA,NA,NA,"list",NA),sep='[.]',remove = FALSE) %>%
  mutate(list = factor(list, 
                       levels=c('any_tiss_hits','cross_hits','protein_hits_filtered','multi_or_cross_hits')),
         list_lab = factor(ifelse(list == 'cross_hits', 'Cross\nExpr',
                                  ifelse(list=='protein_hits_filtered','Cross\nProt',
                                         ifelse(list=='any_tiss_hits','\nWithin\n',
                                                ifelse(list=='multi_or_cross_hits','2+ Within\nor\nWithin + Cross','')))),
                           levels=c('\nWithin\n','Cross\nExpr','Cross\nProt','2+ Within\nor\nWithin + Cross')),
         data_type_lab = factor(ifelse(data_type == 'chip', 'ChIPseq Overlap',
                                       ifelse(data_type == 'motif','Motif Overlap',
                                              ifelse(data_type == 'both','Both Overlap',''))),
                                levels=c('ChIPseq Overlap','Motif Overlap','Both Overlap'))) %>%
  filter(!is.na(list),
         data_type != 'both') %>%
  mutate(p_lab = ifelse(perm_p < 0.0002, "p < 2e-04",
                        ifelse(perm_p < 0.001, paste('p =',format(round(perm_p,4),scientific = FALSE)),
                               ifelse(perm_p < 0.005, paste('p =',format(round(perm_p,3),scientific = FALSE)),
                                      paste('p =',round(perm_p,2))))) ) %>%
  filter(eq=='eq2',
         var == 'true_stat') %>%
  ggplot(aes(list_lab, real)) +
  geom_col(width=.02) +
  geom_hline(yintercept=0) +
  geom_point(col='darkorchid4') +
  geom_point(aes(list_lab,ifelse(real > 0, real*1.5, real)),col=NA) +
  geom_text(aes(list_lab,
                rep(c(0.07,0.0035),4),
                label=p_lab),
            size=3) +
  facet_wrap(~data_type_lab, scales='free_y', nrow=3) +
  theme_classic() +
  theme(strip.text = element_text(hjust=0)) +
  ylab("Overlap Enrichment") +
  xlab('Dataset')
ggsave("plots/fig3_overlap.pdf",
       height=3,width=4)


overlap_sum %>%
  separate(file,into=c(NA,NA,NA,"list",NA),sep='[.]',remove = FALSE) %>%
  mutate(list = factor(list, 
                       levels=c('any_tiss_hits','cross_hits','protein_hits_filtered','multi_or_cross_hits')),
         list_lab = factor(ifelse(list == 'cross_hits', 'Cross\nExpr',
                                  ifelse(list=='protein_hits_filtered','Cross\nProt',
                                         ifelse(list=='any_tiss_hits','\nWithin\n',
                                                ifelse(list=='multi_or_cross_hits','2+ Within\nor\nWithin + Cross','')))),
                           levels=c('\nWithin\n','Cross\nExpr','Cross\nProt','2+ Within\nor\nWithin + Cross')),
         data_type_lab = factor(ifelse(data_type == 'chip', 'ChIPseq Overlap',
                                       ifelse(data_type == 'motif','Motif Overlap',
                                              ifelse(data_type == 'both','Both Overlap',''))),
                                levels=c('ChIPseq Overlap','Motif Overlap','Both Overlap'))) %>%
  filter(!is.na(list),
         data_type != 'both') %>%
  mutate(p_lab = ifelse(perm_p < 0.0002, "p < 2e-04",
                        ifelse(perm_p < 0.001, paste('p =',format(round(perm_p,4),scientific = FALSE)),
                               ifelse(perm_p < 0.005, paste('p =',format(round(perm_p,3),scientific = FALSE)),
                                      paste('p =',round(perm_p,2))))) ) %>%
  filter(eq=='eq2',
         var == 'true_stat') %>%
  ggplot(aes(list_lab, real)) +
  geom_col(width=.02) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=list_lab)) +
  geom_point(aes(list_lab,ifelse(real > 0, real*1.5, real)),col=NA) +
  geom_text(aes(list_lab,
                rep(c(0.07,0.0035),4),
                label=p_lab),
            size=3) +
  facet_wrap(~data_type_lab, scales='free_y', nrow=3) +
  theme_classic() +
  scale_color_manual(values=c('tomato2','gold2','steelblue3','darkorchid3')) +
  theme(strip.text = element_text(hjust=0),
        legend.position = 'none') +
  ylab("Overlap Enrichment") +
  xlab('Dataset')
ggsave("plots/fig3_overlap_col.pdf",
       height=3,width=4)

overlap_sum %>%
  filter(data_type=='chip') %>%
  separate(file,into=c(NA,NA,NA,"list",NA),sep='[.]',remove = FALSE) %>%
  mutate(list = factor(list, 
                       levels=c('any_tiss_hits','cross_hits','protein_hits_filtered','multi_or_cross_hits')),
         list_lab = factor(ifelse(list == 'cross_hits', 'Cross\nExpr',
                                  ifelse(list=='protein_hits_filtered','Cross\nProt',
                                         ifelse(list=='any_tiss_hits','\nWithin\n',
                                                ifelse(list=='multi_or_cross_hits','2+ Within\nor\nWithin + Cross','')))),
                           levels=c('\nWithin\n','Cross\nExpr','Cross\nProt','2+ Within\nor\nWithin + Cross')),
         data_type_lab = factor(ifelse(data_type == 'chip', 'ChIPseq Overlap',
                                       ifelse(data_type == 'motif','Motif Overlap',
                                              ifelse(data_type == 'both','Both Overlap',''))),
                                levels=c('ChIPseq Overlap','Motif Overlap','Both Overlap'))) %>%
  filter(!is.na(list),
         data_type != 'both') %>%
  mutate(p_lab = ifelse(perm_p < 0.0002, "p < 2e-04",
                        ifelse(perm_p < 0.001, paste('p =',format(round(perm_p,4),scientific = FALSE)),
                               ifelse(perm_p < 0.005, paste('p =',format(round(perm_p,3),scientific = FALSE)),
                                      paste('p =',round(perm_p,2))))) ) %>%
  filter(eq=='eq2',
         var == 'true_stat') %>%
  ggplot(aes(list_lab, real)) +
  geom_col(width=.02) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=list_lab)) +
  geom_point(aes(list_lab,ifelse(real > 0, real*1.5, real)),col=NA) +
  geom_text(aes(list_lab,
                rep(c(0.07),4),
                label=p_lab),
            size=3) +
  theme_classic() +
  scale_color_manual(values=c('tomato2','gold2','steelblue3','darkorchid3')) +
  theme(strip.text = element_text(hjust=0),
        legend.position = 'none') +
  ylab("ChIPseq Enrichment") +
  xlab('Dataset')
ggsave("plots/fig3_overlap_col_chip.pdf",
       height=2,width=4)




######## SUPPL FIG STUFF (IRF1)

moc_tf_eq2_chip = read.table("~/projects/overlap/sum_stats3/multi_or_cross_hits.eq2_true.tf_chipseq_stats.txt",
                             header=TRUE,sep='\t')
moc_tf_eq2_mot = read.table("~/projects/overlap/sum_stats3/multi_or_cross_hits.eq2_true.tf_motif_stats.txt",
                             header=TRUE,sep='\t')

merge(moc_tf_eq2_chip, moc_tf_eq2_mot,
      by='tf') %>%
  ggplot(aes(num_lt.x,num_lt.y)) +
  geom_point(aes(col=ifelse(tf%in%c('IRF1','IKZF1'),'red','black'))) +
  geom_text_repel(aes(label=ifelse(tf%in%c('IRF1','IKZF1'),as.character(tf),NA)),
                  nudge_x=4) +
  geom_hline(yintercept=5000) +
  geom_vline(xintercept=5000) +
  theme_classic() +
  xlab('Number of ChIPseq permutations less than real') +
  ylab('Number of motif permutation less than real') +
  theme(legend.position='none') +
  scale_color_manual(values=c('black','forestgreen'))




