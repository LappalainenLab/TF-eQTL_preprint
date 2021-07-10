#!/usr/bin/Rscript
##
##  fig4_plotExpl.R
##
##  EDF 5/12/2021
##

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/projects/MANUSCRIPT/")

temp_data = data.frame(
  cond=factor(c(rep('Promoter',2),rep('Control',2)),levels=c('Promoter','Control')),
  allele=rep(c('ref','alt'),2),
  value=c(36,64,42,36)
)
temp_data %>%
  ggplot(aes(cond,value)) +
  geom_col(aes(alpha=cond,fill=allele),
           position=position_dodge(),
           width=.75) +
  scale_alpha_discrete(range=c(.5,1)) +
  scale_fill_manual(values=c('deepskyblue2','brown3')) +
  theme_classic() +
  ylab('Read Count') +
  xlab('Condition') +
  theme(legend.position='none')
ggsave("plots/fig4_plotExpl1.pdf",
       height=2, width=2)


# temp_data %>%
#   pivot_wider(names_from=allele,values_from=value) %>%
#   mutate(tot_read = ref+alt,
#          alt_af = alt/tot_read) %>%
#   ggplot(aes(cond,alt_af)) +
#   geom_point(aes(color=cond,size=tot_read)) +
#   scale_color_manual(values=c('darkorchid2','darkorchid4')) +
#   scale_size_continuous(limits=c(1,100)) +
#   theme_classic() +
#   ylab('ALT AF') +
#   scale_y_continuous(breaks=c(0,.5,1),
#                      limits=c(0,1)) +
#   xlab('Condition') +
#   theme(legend.position='none')
# ggsave("plots/fig4_plotExpl2.pdf",
#        height=2, width=2)

# 
# mpileups_comb_fi_test_lowp %>%
#   head(n=32) %>%
#   ggplot(aes(1,value)) +
#   geom_col(aes(alpha=cond, fill=allele),
#            position=position_dodge(),
#            width=.75) +
#   geom_text(aes(1,p_height,
#                 label=p_lab)) +
#   scale_alpha_discrete(range=c(.5,1)) +
#   scale_fill_manual(values=c('deepskyblue2','brown3')) +
#   geom_text(aes(1,p_height*1.1,
#                 label=NA)) +
#   facet_wrap(~gene_p_lab,
#              scales = 'free_y') +
#   theme_classic() +
#   ylab('Read Count') +
#   xlab('Time Point')





