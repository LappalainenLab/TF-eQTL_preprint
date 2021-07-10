#!/usr/bin/Rscript
##
##  fig2_plotCrossTiss.R
##
##  EDF 3/29/2021
##


library(dplyr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(ggrepel)
library(ggbeeswarm)

setwd("~/projects/MANUSCRIPT/")

cor_info = read.table("../prot_corrs/tf_measurements/expr_protein_corr.txt",
                       header=TRUE, sep='\t')

expr_hits = read.table("data_tables/cross_tiss_corrs/cross_corrs_fdr05.txt",
                       header=TRUE, sep='\t')
prot_hits = read.table("../overlap/replication/sig_assoc.fdr20.protein_hits_filtered.txt",
                       header=TRUE, sep='\t')

prot_hits_sum = prot_hits %>%
  group_by(tf) %>%
  summarise(num_sig = n(),
            data='prot')

expr_hits_sum = expr_hits %>%
  group_by(tf) %>%
  summarise(num_sig = n(),
            data='expr')

both_hits_sum = rbind(prot_hits_sum, expr_hits_sum) %>%
  merge(cor_info, by.x='tf', by.y='description')

both_hits_sum_plotting = both_hits_sum %>%
  pivot_wider(id_cols=c('tf','sp_rho'),
              names_from='data',
              values_from='num_sig') %>%
  mutate(prot=ifelse(is.na(prot),-500,prot),
         sp_rho_plot=ifelse(prot==-500,NA,sp_rho))

both_hits_sum_plotting %>%
  ggplot(aes(prot,expr)) +
  geom_text_repel(aes(label=ifelse(expr>8000|prot>1000,as.character(tf),'')),
                  size=3,nudge_x = 200) +
  geom_point(aes(fill=sp_rho_plot),shape=21,size=2.5) +
  theme_classic() +
  ylab("Expression-based at 5% FDR") +
  xlab("Protein-based at 20% FDR") +
  theme(legend.position = c(0.8, 0.7),
        legend.text = element_text(size=6),
        legend.title = element_text(size=9)) +
  scale_fill_gradient2(low='blue',mid='white',high='red', 
                       guide = guide_colorbar(barwidth = .5, barheight = 3)) +
  labs(fill=element_text("Expr:Prot rho")) +
  scale_x_continuous(breaks=c(-500,0,1000,2000,3000),
                     labels=c("NA", "0", "1000",'2000','3000'),
                     limits = c(-500,3800))
ggsave("plots/fig2_cross_tiss.pdf",
       height=3, width=3.5)

# 
# both_hits_sum_plotting %>%
#   ggplot(aes(expr,prot)) +
#   geom_text_repel(aes(label=ifelse(expr>8000|prot>1000,as.character(tf),'')),
#                   size=3,nudge_x = 200) +
#   geom_point(aes(fill=sp_rho_plot),shape=21,size=2.5) +
#   theme_classic() +
#   xlab("Expression-based at 5% FDR") +
#   ylab("Protein-based at 20% FDR") +
#   theme(legend.position = c(0.8, 0.7),
#         legend.text = element_text(size=6),
#         legend.title = element_text(size=9)) +
#   scale_fill_gradient2(low='blue',mid='white',high='red', 
#                        guide = guide_colorbar(barwidth = .5, barheight = 3)) +
#   labs(fill=element_text("Expr:Prot rho")) +
#   scale_y_continuous(breaks=c(-500,0,1000,2000,3000),
#                      labels=c("NA", "0", "1000",'2000','3000'),
#                      limits = c(-500,3800))


both_hits_sum %>%
  mutate(xnames=ifelse(data=='expr',"Expression\n5% FDR",
                       "Protein\n20% FDR")) %>%
  ggplot(aes(xnames,num_sig)) +
  geom_violin(aes(y=ifelse(data=='prot'&num_sig>2000,NA,num_sig))) +
  geom_point(aes(y=ifelse(data=='expr'&num_sig<=7000,num_sig,NA)),
             position=position_jitter(width=.2),
             size=.5) +
  geom_point(aes(y=ifelse(data=='expr'&num_sig>7000,num_sig,NA)),
             position=position_jitter(width=.02),
             size=.5) +
  geom_point(aes(y=ifelse(data=='prot'&num_sig<1000,num_sig,NA)),
             position=position_jitter(width=.2),
             size=.5) +
  geom_point(aes(y=ifelse(data=='prot'&num_sig>=1000,num_sig,NA)),
             position=position_jitter(width=0),
             size=.5) +
  geom_text_repel(aes(label=ifelse(data=='expr'&num_sig>8000,as.character(tf),
                                   ifelse(data=='prot'&num_sig>1000,as.character(tf),
                                          '')))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12)) +
  ylab('Cross-tissue TF-eQTLs')
both_hits_sum %>%
  mutate(xnames=ifelse(data=='expr',"Expression\n5% FDR",
                       "Protein\n20% FDR")) %>%
  ggplot(aes(xnames,num_sig)) +
  geom_quasirandom(aes(xnames,num_sig),
             size=.5) +
  geom_text_repel(aes(label=ifelse(data=='expr'&num_sig>8000,as.character(tf),
                                   ifelse(data=='prot'&num_sig>1000,as.character(tf),
                                          '')))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12)) +
  ylab('Cross-tissue TF-eQTLs')
ggsave("plots/fig2_cross_tiss2.pdf",
       height=3, width=3)
