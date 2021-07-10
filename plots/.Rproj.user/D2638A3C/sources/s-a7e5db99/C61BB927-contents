#!/usr/bin/Rscript
##
##  fig2_plotWinTiss.R
##
##  EDF 3/29/2021
##

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(forcats)


setwd("~/projects/MANUSCRIPT/")

tf_info = read.table("../TFi-eQTL/input_files/TFs.info.curated.txt",
                     header=TRUE, sep='\t')
tiss_info = read.table("../TFi-eQTL/input_files/tissue_info.final.txt",
                       header=TRUE, sep='\t') %>%
  mutate(tiss_type = factor(
    c('fat','epithelial','brain','brain','brain',
      'blood','muscle','organ','epithelial','muscle',
      'nervous','organ','nervous','epithelial','blood','nervous'),
    levels=c('blood','fat','brain','nervous','epithelial','muscle','organ')
  ),
  tiss_type_2 = factor(ifelse(tiss_type=='epithelial','epithelial',
                       ifelse(tiss_type=='nervous','nervous',
                              'internal'))))

# sig_hits = read.table("data_tables/within_tiss_corrs/win_corrs_fdr20.txt",
#                       header=TRUE, sep='\t')
sig_hits = read.table("data_tables/within_tiss_corrs/win_corrs_fdr20.sig_eqtl.txt",
                      header=TRUE, sep='\t') %>%
  merge(tf_info %>% select(1:2) %>% rename(tf_id = Name),
        by.x=c('tf'), by.y=c('TF')) %>%
  merge(tiss_info, by.x="tiss", by.y="TISSUE_NAME") %>%
  filter(tf_id != as.character(gene)) %>%
  mutate(tiss = factor(tiss, levels=unique(as.character(tiss_info$TISSUE_NAME))))

top_per_tiss = sig_hits %>%
  group_by(TISSUE_ABBRV,tf) %>%
  summarize(count=n()) %>%
  group_by(TISSUE_ABBRV) %>%
  summarize(top_tf = tf[which.max(count)])

sig_hits %>%
  ggplot(aes(TISSUE_ABBRV)) +
  geom_bar(aes(fill=tiss)) +
  theme_classic() +
  scale_fill_manual(values = as.character(tiss_info$TISSUE_RCOL)) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        axis.text.y = element_text(size=12)) +
  xlab('Tissue') +
  ylab('Within-tissue at 20% FDR') +
  geom_text(data=top_per_tiss, 
            aes(TISSUE_ABBRV,rep(c(7000,6000),8),
                label=top_tf),
            size=3) +
  scale_y_continuous(breaks=c(0,1000,2000,3000,4000,5000, 6500),
                     labels=c("0", "1000", "2000",'3000','4000','5000','Top TF'))
ggsave("plots/fig2_win_tiss.pdf",
       height=3,width=6)



sig_by_tf = sig_hits %>%
  group_by(TISSUE_ABBRV,tiss,tf) %>%
  summarize(count=n()) %>%
  group_by(TISSUE_ABBRV) %>%
  mutate(max_tf = tf[which.max(count)])

sig_by_tf %>%
  ggplot(aes(TISSUE_ABBRV,count)) +
  geom_violin() +
  geom_point(aes(col=tiss),
             position=position_jitter(width=.1),
             size=.5) +
  scale_color_manual(values = as.character(tiss_info$TISSUE_RCOL)) +
  ylim(0,210) +
  geom_text_repel(aes(label=ifelse(tf==max_tf,as.character(tf),'')),
                  nudge_y = 10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        axis.text.y = element_text(size=12)) +
  xlab('Tissue') +
  ylab('Within-tissue at 20% FDR')
ggsave("plots/fig2_win_tiss2.pdf",
       height=3,width=6)

sig_by_tf %>%
  merge(tiss_info) %>%
  ggplot(aes(fct_reorder(TISSUE_ABBRV,as.numeric(tiss_type)),count)) +
  geom_violin() +
  geom_point(aes(col=tiss_type),
             position=position_jitter(width=.1),
             size=.5) +
  ylim(0,210) +
  geom_text_repel(aes(label=ifelse(tf==max_tf,as.character(tf),'')),
                  nudge_y = 10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        axis.text.y = element_text(size=12)) +
  xlab('Tissue') +
  ylab('Within-tissue at 20% FDR') +
  scale_color_manual(values=c('hotpink2','chocolate2','goldenrod2',
                              'olivedrab3','deepskyblue3','slateblue3',
                              'antiquewhite3'))
ggsave("plots/fig2_win_tiss2_orgcol.pdf",
       height=3,width=6)
