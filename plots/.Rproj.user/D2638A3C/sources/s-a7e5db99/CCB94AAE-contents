#!/usr/bin/Rscript
##
##  fig1_plotTheory.R
##  
##  EDF 5/4/2021
##


library(ggplot2)
library(reshape2)
library(dplyr)

## generate theoretical plot

sigmoid = function(x) { 1/ (1+exp(-x)) }
sigmoid_2 = function(x) { 1/ (1+.1*exp(-x)) }
x <- seq(-8,10,.01)
data.frame(x,sigmoid(x),sigmoid_2(x)) -> sigmoid_df
names(sigmoid_df) <- c('x','sig1','sig2')

ggplot(sigmoid_df,aes(x)) +
  geom_line(aes(y=sig1),color='skyblue3',size=2) +
  geom_line(aes(y=sig2),color='darkred',size=2) + 
  theme_classic() +
  xlab(label="TF level") +
  ylab(label="TF binding occupancy / 
target gene expression") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
      axis.title = element_text(size=14))
ggsave("fig1_theory1.pdf",
       width=3, height=3)

# names(sigmoid_df) <- c('x','alt allele','ref allele')
# melt(sigmoid_df,id="x") -> sigmoid_df_long
# names(sigmoid_df_long) <- c('x','genotype','sig')
# ggplot(sigmoid_df_long,aes(x,sig)) +
#   geom_line(size=2,aes(color=genotype)) + 
#   theme_classic() +
#   xlab(label="TF level") +
#   ylab(label="TF binding occupancy / 
# target gene expression") +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_text(size=20),
#         legend.position=c(.8,.4),
#         legend.background=element_rect(linetype="solid",color="black"),
#         legend.title=element_blank(),
#         legend.text=element_text(size=20)) +
#   scale_color_manual(values=c('red','black'))

low=100
high=10000
bg=.5
kd = c(low,high)
TF_expr = seq(0,6,.1)
occ_by_kd = as.data.frame(sapply(kd, function(kd) { sapply(10^TF_expr, function(x) { (x/(kd+x))*(1-bg)+bg })}))
names(occ_by_kd) <- kd
occ_by_kd$TF_expr <- 10^TF_expr
occ_by_kd_long = melt(occ_by_kd, id=c('TF_expr'))
names(occ_by_kd_long) <- c('TF_expr','Kd','occ')

ggplot(occ_by_kd_long,aes(log(TF_expr),occ)) +
  geom_line(size=2,aes(color=Kd)) + 
  theme_classic() +
  xlab(label="log(TF level)") +
  ylab(label="DNA occupancy / 
target gene expression") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=20),
        legend.position=c(.8,.4),
        legend.background=element_rect(linetype="solid",color="black"),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20, hjust=0.5))

occ_by_kd$aFC <- unlist(log2(occ_by_kd[c(as.character(low))]/occ_by_kd[c(as.character(high))]))
ggplot(occ_by_kd,aes(log(TF_expr),aFC)) +
  geom_line(size=2) + 
  theme_classic() +
  xlab(label="log(TF level)") +
  ylab(label="eQTL effect size") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size=14))
ggsave("fig1_theory2.pdf",
       width=3, height=3)

sapply(c(10^-25,.00001,.0001,.001,.01,.1,.5,1,2,3,4,5,6,7,8,9,10,20,30,40,50),function(x) { (x/(1+x))+.01 }) -> ys_1kd


set.seed(9342)
data.frame(x=seq(1:15)+rnorm(15,sd = .1),
           y=seq(1:15)/5+rnorm(15)) %>%
  ggplot(aes(x,y)) +
  geom_point()  + 
  geom_smooth(method='lm', 
              formula= y~x,
              se=FALSE,
              col='black') +
  theme_classic() +
  xlab(label="TF level") +
  ylab(label="eQTL effect size") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size=14))
ggsave("fig1_theory3.pdf",
       width=3, height=3)





