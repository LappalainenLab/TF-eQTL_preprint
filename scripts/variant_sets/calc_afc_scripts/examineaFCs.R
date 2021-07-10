#!/usr/bin/Rscript
##
##  examineaFCs.R
##
##  EDF 4/23/2020
##

library(reshape2)
library(ggplot2)

setwd("~/projects/crosstiss_tf_corrs/afcs/")

tiss_names=read.table("../input_files/tissue_translation_colors_v8.txt",
                      header=TRUE, sep='\t')
tiss_counts=read.table("../input_files/tissue_sample_count.txt",
                       header=FALSE,sep='\t')
names(tiss_counts) <- c('TISSUE_NAME','num_samps')
tiss_trans=merge(tiss_names, tiss_counts,
                 by=c('TISSUE_NAME'))

afcs_cov = read.table("combined/all_tiss.afcs.cov.overlap_vars.MAF05.txt.gz",
                      header=TRUE, sep = '\t')
afcs_cov_long = melt(afcs_cov, id.vars=c('sid','chr','pos','pid'),
                     variable='tiss', value.name='afc')

# ggplot(afcs_cov_long) +
#   facet_wrap(vars(tiss)) +
#   geom_histogram(aes(afc)) +
#   theme_classic
# 
# ggplot(afcs_cov_long) +
#   geom_violin(aes(tiss, afc)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) 

afcs_cov_long %>%
  mutate(tiss = factor(tiss, 
    levels=as.character(tiss_trans$TISSUE_ABBRV[
    order(-tiss_trans$num_samps)])) ) %>%
  ggplot(aes(tiss, afc)) +
  geom_violin(aes(fill=tiss)) +
  scale_fill_manual(values=as.character(tiss_trans$TISSUE_RCOL[
    order(-tiss_trans$num_samps)])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  theme(legend.position = "none")

afcs_cov_tiss_sum = afcs_cov_long %>%
  group_by(tiss) %>%
  summarize(mean=mean(afc, na.rm=TRUE), 
            stdev=sd(afc, na.rm=TRUE)) %>%
  merge(tiss_trans[c(3,4,6)],
        by.x='tiss',by.y='TISSUE_ABBRV') %>%
  arrange(tiss)
afcs_cov_tiss_sum %>%
  ggplot(aes(num_samps,stdev)) +
  geom_point(aes(col=tiss)) +
  scale_color_manual(values=as.character(afcs_cov_tiss_sum$TISSUE_RCOL)) +
  theme_classic() +
  xlim(c(0,800)) +
  ylim(c(0.3,0.62)) +
  theme(legend.position = 'none')



afcs_nocov = read.table("combined/all_tiss.afcs.nocov.overlap_vars.MAF05.txt.gz",
                      header=TRUE, sep = '\t')
afcs_nocov_long = melt(afcs_nocov, id.vars=c('sid','chr','pos','pid'),
                     variable='tiss', value.name='afc')
afcs_nocov_long = afcs_nocov_long %>%
  filter(!(tiss %in% c("CVXECT","CVSEND","FLLPNT","KDNMDL")) &
           !(is.na(tiss)))

# ggplot(afcs_nocov_long) +
#   facet_wrap(vars(tiss)) +
#   geom_histogram(aes(afc)) +
#   theme_classic
# 
# ggplot(afcs_nocov_long) +
#   geom_violin(aes(tiss, afc)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) 

afcs_nocov_long %>%
  mutate(tiss = factor(tiss, 
                       levels=as.character(tiss_trans$TISSUE_ABBRV[
                         order(-tiss_trans$num_samps)])) ) %>%
  ggplot(aes(tiss, afc)) +
  geom_violin(aes(fill=tiss)) +
  scale_fill_manual(values=as.character(tiss_trans$TISSUE_RCOL[
    order(-tiss_trans$num_samps)])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  theme(legend.position = "none")

afcs_nocov_tiss_sum = afcs_nocov_long %>%
  group_by(tiss) %>%
  summarize(mean=mean(afc, na.rm=TRUE), 
            stdev=sd(afc, na.rm=TRUE)) %>%
  merge(tiss_trans[c(3,4,6)],
        by.x='tiss',by.y='TISSUE_ABBRV') %>%
  arrange(tiss)
afcs_nocov_tiss_sum %>%
  ggplot(aes(num_samps,stdev)) +
  geom_point(aes(col=tiss)) +
  scale_color_manual(values=as.character(afcs_nocov_tiss_sum$TISSUE_RCOL)) +
  theme_classic() +
  xlim(c(0,800)) +
  ylim(c(0.3,0.62)) +
  theme(legend.position = 'none')
