#!/usr/bin/Rscript
##
##  fig4_mocexamps.R
##
##  EDF 5/12/2021
##

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/projects/MANUSCRIPT/")

fi_tst_moc = read.table("data_tables/IRF1_kd/Fi_tst_combined_byvar_timecomb.mocgenes.txt",
                        header=TRUE,sep='\t')

fi_tst_moc %>%
  pivot_longer(cols=c(altp_C,altp_P,tot_C,tot_P)) %>%
  separate(name,sep='_',into=c('stat','cond'),) %>%
  pivot_wider(names_from=stat) %>%
  mutate(cond=factor(cond,levels=c('P','C'))) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = .95) %>%
  group_by(gene,description,exon,variant,chr,pos,ref,alt) %>%
  mutate(gene_p_lab = paste(description,"p =",
                            ifelse(fi_p<0.01,round(fi_p,3),round(fi_p,2)))) %>%
  ggplot(aes(1,altp)) +
  geom_point(aes(alpha=cond,group=cond,
                 size=tot),
           position=position_dodge(.75),
           color='forestgreen') +
  geom_text(aes(1,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  facet_wrap(~gene_p_lab,
             scales = 'free_y',
             nrow = 2) +
  theme_classic() +
  ylim(0,1) +
  ylab('Read Count') +
  xlab('Time Point')




fi_tst_moc %>%
  pivot_longer(cols=c(altp_C,altp_P,tot_C,tot_P)) %>%
  separate(name,sep='_',into=c('stat','cond'),) %>%
  pivot_wider(names_from=stat) %>%
  mutate(cond=factor(cond,levels=c('P','C'))) %>%
  mutate(p_lab = ifelse(fi_p < 0.01, "**",
                        ifelse(fi_p < 0.05, "*",
                               ifelse(fi_p < 0.1, "+",""))),
         p_height = .95) %>%
  ggplot() +
  geom_segment(data=fi_tst_moc,aes(x=description,xend=description,
                                   y=altp_P,yend=altp_C)) +
  geom_point(aes(description,altp,
                 alpha=cond,group=cond,
                 size=tot),
             color='forestgreen') +
  geom_text(aes(description,p_height,
                label=p_lab)) +
  scale_alpha_discrete(range=c(.5,1)) +
  theme_classic() +
  ylim(0,1) +
  ylab('ALT AF') +
  xlab('Gene') +
  coord_flip()
ggsave("plots/fig4_mocexamps.pdf",
       width=4,height=2)


