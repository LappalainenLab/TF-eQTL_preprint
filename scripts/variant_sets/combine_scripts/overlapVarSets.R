#!/usr/bin/Rscript
##
##  overlapVarSets.R
##
##  EDF 4/6/2020
##

library(VennDiagram)
library(RColorBrewer)
library(dplyr)

setwd("~/projects/TFi-eQTL/")


cav_vars_05 = read.table("variant_sets/caviar/caviar_var.95set.all_tiss.MAF05.vcf.gz", 
                         header=FALSE, 
                         colClasses = c("character", "integer","character", rep("NULL", 5)))
cav_vars_10 = read.table("variant_sets/caviar/caviar_var.95set.all_tiss.MAF10.vcf.gz",
                         header=FALSE, 
                         colClasses = c("character", "integer","character", rep("NULL", 5)))

mot_either_05_vars = read.table("variant_sets/HOCOMOCO/HOCOMOCO_TF_motif_overlap.either.GTEx_v8.curated_set.MAF05.vars.list",
                                colClasses = c('character'))
mot_change_05_vars = read.table("variant_sets/HOCOMOCO/HOCOMOCO_TF_motif_overlap.change.GTEx_v8.curated_set.MAF05.vars.list",
                                colClasses = c('character'))
mot_either_10_vars = read.table("variant_sets/HOCOMOCO/HOCOMOCO_TF_motif_overlap.either.GTEx_v8.curated_set.MAF10.vars.list",
                                colClasses = c('character'))
mot_change_10_vars = read.table("variant_sets/HOCOMOCO/HOCOMOCO_TF_motif_overlap.change.GTEx_v8.curated_set.MAF10.vars.list",
                                colClasses = c('character'))

chip_vars_05 = read.table("variant_sets/ENCODE/ENCODE_TF_ChIPseq_overlap.GTEx_v8.curated_set.MAF05.vars.list",
                          colClasses = c('character'))
chip_vars_10 = read.table("variant_sets/ENCODE/ENCODE_TF_ChIPseq_overlap.GTEx_v8.curated_set.MAF10.vars.list",
                          colClasses = c('character'))


## 5% maf overlap (with either mot)
all_05 = 6539591
cols = brewer.pal(3, "Pastel2")
venn.diagram( x = list( cav_vars_05$V3, 
                       mot_either_05_vars$V1, 
                       chip_vars_05$V1 ),
              category.names = c( paste0('caviar 95% set\nn=',nrow(cav_vars_05)),
                                 paste0('motif in either\nn=',nrow(mot_either_05_vars)),
                                 paste0('chipseq\nn=',nrow(chip_vars_05)) ),
              filename="variant_sets/plots/overlap_05.png",
              imagetype='png',
              total.population = all_05,
              lty = 'blank',
              fill = cols,
              cex = 1,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.fontfamily = "sans",
              cat.pos = c(-27, 27, 135),
              width=2400,
              height=2400)
both = sum(cav_vars_05$V3 %in% mot_either_05_vars$V1)
fisher.test( matrix( c( both,
                     nrow(cav_vars_05)-both,
                     nrow(mot_either_05_vars)-both,
                     all_05 - nrow(cav_vars_05) - nrow(mot_either_05_vars) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(mot_either_05_vars)-both,
                        all_05 - nrow(cav_vars_05) - nrow(mot_either_05_vars) + both),
                     nrow =2 ) )$p.value
both = sum(cav_vars_05$V3 %in% chip_vars_05$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(cav_vars_05) - nrow(chip_vars_05) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(cav_vars_05) - nrow(chip_vars_05) + both),
                     nrow =2 ) )$p.value
both = sum(mot_either_05_vars$V1 %in% chip_vars_05$V1)
fisher.test( matrix( c( both,
                        nrow(mot_either_05_vars)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(mot_either_05_vars) - nrow(chip_vars_05) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(mot_either_05_vars)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(mot_either_05_vars) - nrow(chip_vars_05) + both),
                     nrow =2 ) )$p.value

all = sum(mot_either_05_vars$V1 %in% chip_vars_05$V1 &
            mot_either_05_vars$V1 %in% cav_vars_05$V3)
fisher.test( matrix( c( all,
                        sum(chip_vars_05$V1 %in% cav_vars_05$V3)-all,
                        nrow(mot_either_05_vars)-all,
                        all_05 - sum(chip_vars_05$V1 %in% cav_vars_05$V3) - nrow(mot_either_05_vars) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(chip_vars_05$V1 %in% cav_vars_05$V3)-all,
                        nrow(mot_either_05_vars)-all,
                        all_05 - sum(chip_vars_05$V1 %in% cav_vars_05$V3) - nrow(mot_either_05_vars) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_either_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_either_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_either_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_either_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_either_05_vars$V1 %in% chip_vars_05$V1)-all,
                        nrow(cav_vars_05)-all,
                        all_05 - sum(mot_either_05_vars$V1 %in% chip_vars_05$V1) - nrow(cav_vars_05) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_either_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_either_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )$p.value

# write.table(mot_either_05_vars[mot_either_05_vars$V1 %in% chip_vars_05$V1 &
#                                  mot_either_05_vars$V1 %in% cav_vars_05$V3, ], 
#             "variant_sets/overlap_vars.MAF05.list",
#             row.names=FALSE, quote=FALSE, col.names=FALSE)


## 05% maf overlap (with change motif)
all_05 = 6539591
cols = brewer.pal(3, "Pastel2")
venn.diagram( x = list( cav_vars_05$V3, 
                        mot_change_05_vars$V1, 
                        chip_vars_05$V1 ),
              category.names = c( paste0('caviar 95% set\nn=',nrow(cav_vars_05)),
                                  paste0('motif change\nn=',nrow(mot_change_05_vars)),
                                  paste0('chipseq\nn=',nrow(chip_vars_05)) ),
              filename="variant_sets/plots/overlap_05_change.png",
              imagetype='png',
              total.population = all_05,
              lty = 'blank',
              fill = cols,
              cex = 1,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.fontfamily = "sans",
              cat.pos = c(-27, 27, 135),
              width=2400,
              height=2400)
both = sum(cav_vars_05$V3 %in% mot_change_05_vars$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(mot_change_05_vars)-both,
                        all_05 - nrow(cav_vars_05) - nrow(mot_change_05_vars) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(mot_change_05_vars)-both,
                        all_05 - nrow(cav_vars_05) - nrow(mot_change_05_vars) + both),
                     nrow =2 ) )$p.value
both = sum(cav_vars_05$V3 %in% chip_vars_05$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(cav_vars_05) - nrow(chip_vars_05) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_05)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(cav_vars_05) - nrow(chip_vars_05) + both),
                     nrow =2 ) )$p.value
both = sum(mot_change_05_vars$V1 %in% chip_vars_05$V1)
fisher.test( matrix( c( both,
                        nrow(mot_change_05_vars)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(mot_change_05_vars) - nrow(chip_vars_05) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(mot_change_05_vars)-both,
                        nrow(chip_vars_05)-both,
                        all_05 - nrow(mot_change_05_vars) - nrow(chip_vars_05) + both),
                     nrow =2 ) )$p.value

all = sum(mot_change_05_vars$V1 %in% chip_vars_05$V1 &
            mot_change_05_vars$V1 %in% cav_vars_05$V3)
fisher.test( matrix( c( all,
                        sum(chip_vars_05$V1 %in% cav_vars_05$V3)-all,
                        nrow(mot_change_05_vars)-all,
                        all_05 - sum(chip_vars_05$V1 %in% cav_vars_05$V3) - nrow(mot_change_05_vars) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(chip_vars_05$V1 %in% cav_vars_05$V3)-all,
                        nrow(mot_change_05_vars)-all,
                        all_05 - sum(chip_vars_05$V1 %in% cav_vars_05$V3) - nrow(mot_change_05_vars) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_change_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_change_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_change_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_change_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_change_05_vars$V1 %in% chip_vars_05$V1)-all,
                        nrow(cav_vars_05)-all,
                        all_05 - sum(mot_change_05_vars$V1 %in% chip_vars_05$V1) - nrow(cav_vars_05) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_change_05_vars$V1 %in% cav_vars_05$V3)-all,
                        nrow(chip_vars_05)-all,
                        all_05 - sum(mot_change_05_vars$V1 %in% cav_vars_05$V3) - nrow(chip_vars_05) + all),
                     nrow =2 ) )$p.value

## 10% maf overlaps (with either motif)
all_10 = 5016647
cols = brewer.pal(3, "Pastel2")
venn.diagram( x = list( cav_vars_10$V3, 
                        mot_either_10_vars$V1, 
                        chip_vars_10$V1 ),
              category.names = c( paste0('caviar 95% set\nn=',nrow(cav_vars_10)),
                                  paste0('motif in either\nn=',nrow(mot_either_10_vars)),
                                  paste0('chipseq\nn=',nrow(chip_vars_10)) ),
              filename="variant_sets/plots/overlap_10.png",
              imagetype='png',
              total.population = all_10,
              lty = 'blank',
              fill = cols,
              cex = 1,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.fontfamily = "sans",
              cat.pos = c(-27, 27, 135),
              width=2400,
              height=2400)
both = sum(cav_vars_10$V3 %in% mot_either_10_vars$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(mot_either_10_vars)-both,
                        all_10 - nrow(cav_vars_10) - nrow(mot_either_10_vars) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(mot_either_10_vars)-both,
                        all_10 - nrow(cav_vars_10) - nrow(mot_either_10_vars) + both),
                     nrow =2 ) )$p.value
both = sum(cav_vars_10$V3 %in% chip_vars_10$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(cav_vars_10) - nrow(chip_vars_10) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(cav_vars_10) - nrow(chip_vars_10) + both),
                     nrow =2 ) )$p.value
both = sum(mot_either_10_vars$V1 %in% chip_vars_10$V1)
fisher.test( matrix( c( both,
                        nrow(mot_either_10_vars)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(mot_either_10_vars) - nrow(chip_vars_10) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(mot_either_10_vars)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(mot_either_10_vars) - nrow(chip_vars_10) + both),
                     nrow =2 ) )$p.value

all = sum(mot_either_10_vars$V1 %in% chip_vars_10$V1 &
            mot_either_10_vars$V1 %in% cav_vars_10$V3)
fisher.test( matrix( c( all,
                        sum(chip_vars_10$V1 %in% cav_vars_10$V3)-all,
                        nrow(mot_either_10_vars)-all,
                        all_10 - sum(chip_vars_10$V1 %in% cav_vars_10$V3) - nrow(mot_either_10_vars) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(chip_vars_10$V1 %in% cav_vars_10$V3)-all,
                        nrow(mot_either_10_vars)-all,
                        all_10 - sum(chip_vars_10$V1 %in% cav_vars_10$V3) - nrow(mot_either_10_vars) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_either_10_vars$V1 %in% cav_vars_10$V3)-all,
                        nrow(chip_vars_10)-all,
                        all_10 - sum(mot_either_10_vars$V1 %in% cav_vars_10$V3) - nrow(chip_vars_10) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_either_10_vars$V1 %in% cav_vars_10$V3)-all,
                        nrow(chip_vars_10)-all,
                        all_10 - sum(mot_either_10_vars$V1 %in% cav_vars_10$V3) - nrow(chip_vars_10) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_either_10_vars$V1 %in% chip_vars_10$V1)-all,
                        nrow(cav_vars_10)-all,
                        all_10 - sum(mot_either_10_vars$V1 %in% chip_vars_10$V1) - nrow(cav_vars_10) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_either_10_vars$V1 %in% chip_vars_10$V1)-all,
                        nrow(cav_vars_10)-all,
                        all_10 - sum(mot_either_10_vars$V1 %in% chip_vars_10$V1) - nrow(cav_vars_10) + all),
                     nrow =2 ) )$p.value
# write.table(mot_either_10_vars[mot_either_10_vars$V1 %in% chip_vars_10$V1 &
#                                  mot_either_10_vars$V1 %in% cav_vars_10$V3, ], 
#             "variant_sets/overlap_vars.MAF10.list",
#             row.names=FALSE, quote=FALSE, col.names=FALSE)



## 10% maf overlap (with change motif)
all_10 = 5016647
cols = brewer.pal(3, "Pastel2")
venn.diagram( x = list( cav_vars_10$V3, 
                        mot_change_10_vars$V1, 
                        chip_vars_10$V1 ),
              category.names = c( paste0('caviar 95% set\nn=',nrow(cav_vars_10)),
                                  paste0('motif change\nn=',nrow(mot_change_10_vars)),
                                  paste0('chipseq\nn=',nrow(chip_vars_10)) ),
              filename="variant_sets/plots/overlap_10_change.png",
              imagetype='png',
              total.population = all_10,
              lty = 'blank',
              fill = cols,
              cex = 1,
              fontface = "bold",
              fontfamily = "sans",
              cat.cex = 1,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.fontfamily = "sans",
              cat.pos = c(-27, 27, 135),
              width=2400,
              height=2400)
both = sum(cav_vars_10$V3 %in% mot_change_10_vars$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(mot_change_10_vars)-both,
                        all_10 - nrow(cav_vars_10) - nrow(mot_change_10_vars) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(mot_change_10_vars)-both,
                        all_10 - nrow(cav_vars_10) - nrow(mot_change_10_vars) + both),
                     nrow =2 ) )$p.value
both = sum(cav_vars_10$V3 %in% chip_vars_10$V1)
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(cav_vars_10) - nrow(chip_vars_10) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(cav_vars_10)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(cav_vars_10) - nrow(chip_vars_10) + both),
                     nrow =2 ) )$p.value
both = sum(mot_change_10_vars$V1 %in% chip_vars_10$V1)
fisher.test( matrix( c( both,
                        nrow(mot_change_10_vars)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(mot_change_10_vars) - nrow(chip_vars_10) + both),
                     nrow =2 ) )
fisher.test( matrix( c( both,
                        nrow(mot_change_10_vars)-both,
                        nrow(chip_vars_10)-both,
                        all_10 - nrow(mot_change_10_vars) - nrow(chip_vars_10) + both),
                     nrow =2 ) )$p.value

all = sum(mot_change_10_vars$V1 %in% chip_vars_10$V1 &
            mot_change_10_vars$V1 %in% cav_vars_10$V3)
fisher.test( matrix( c( all,
                        sum(chip_vars_10$V1 %in% cav_vars_10$V3)-all,
                        nrow(mot_change_10_vars)-all,
                        all_10 - sum(chip_vars_10$V1 %in% cav_vars_10$V3) - nrow(mot_change_10_vars) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(chip_vars_10$V1 %in% cav_vars_10$V3)-all,
                        nrow(mot_change_10_vars)-all,
                        all_10 - sum(chip_vars_10$V1 %in% cav_vars_10$V3) - nrow(mot_change_10_vars) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_change_10_vars$V1 %in% cav_vars_10$V3)-all,
                        nrow(chip_vars_10)-all,
                        all_10 - sum(mot_change_10_vars$V1 %in% cav_vars_10$V3) - nrow(chip_vars_10) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_change_10_vars$V1 %in% cav_vars_10$V3)-all,
                        nrow(chip_vars_10)-all,
                        all_10 - sum(mot_change_10_vars$V1 %in% cav_vars_10$V3) - nrow(chip_vars_10) + all),
                     nrow =2 ) )$p.value

fisher.test( matrix( c( all,
                        sum(mot_change_10_vars$V1 %in% chip_vars_10$V1)-all,
                        nrow(cav_vars_10)-all,
                        all_10 - sum(mot_change_10_vars$V1 %in% chip_vars_10$V1) - nrow(cav_vars_10) + all),
                     nrow =2 ) )
fisher.test( matrix( c( all,
                        sum(mot_change_10_vars$V1 %in% chip_vars_10$V1)-all,
                        nrow(cav_vars_10)-all,
                        all_10 - sum(mot_change_10_vars$V1 %in% chip_vars_10$V1) - nrow(cav_vars_10) + all),
                     nrow =2 ) )$p.value


## get eqtls (genes) from caviar list
caviar_eqtls=read.table("variant_sets/caviar/caviar_var.uniqueeqtls.95set.txt",
                        header=FALSE, sep='\t')
#caviar_eqtls$chr_pos = paste(caviar_eqtls$V2, caviar_eqtls$V3,sep='_')
cav_overlap_05 = cav_vars_05[cav_vars_05$V3 %in% mot_either_05_vars$V1 &
                               cav_vars_05$V3 %in% chip_vars_05$V1, ]
#cav_overlap_05$chr_pos = paste(cav_overlap_05$V1, cav_overlap_05$V2, sep='_')
cav_overlap_10 = cav_vars_10[cav_vars_10$V3 %in% mot_either_10_vars$V1 &
                               cav_vars_10$V3 %in% chip_vars_10$V1, ]
#cav_overlap_10$chr_pos = paste(cav_overlap_10$V1, cav_overlap_10$V2, sep='_')

cav_eqtls_05 = merge(caviar_eqtls, cav_overlap_05,
                     by.x=c('V2','V3'), by.y=c('V1','V2'))
names(cav_eqtls_05) <- c('chr','pos','gene','var')
cav_eqtls_05$chr <- factor(cav_eqtls_05$chr,
                              levels=paste0('chr',seq(1,22)))
cav_eqtls_05$pos <- as.integer(cav_eqtls_05$pos)
cav_eqtls_05_gene_sorted = cav_eqtls_05 %>%
  arrange(gene,chr,pos) %>%
  select(gene, chr, pos, var)
write.table(cav_eqtls_05_gene_sorted, "variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
            col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
sapply(paste0('chr',seq(1,22)), function(chri) { 
  cav_eqtls_05_gene_sorted %>% filter(chr==chri) %>%
    write.table(paste0("variant_sets/by_chr/caviar_var.95set.eqtls.MAF05.overlap.genesort.",chri,".txt"),
                col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
})
cav_eqtls_05_var_sorted = cav_eqtls_05 %>%
  arrange(chr,pos, gene) %>%
  select(chr, pos, gene, var)
write.table(cav_eqtls_05_var_sorted, "variant_sets/caviar_var.95set.eqtls.MAF05.overlap.varsort.txt",
            col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
sapply(paste0('chr',seq(1,22)), function(chri) { 
  cav_eqtls_05_var_sorted %>% filter(chr==chri) %>%
    write.table(paste0("variant_sets/by_chr/caviar_var.95set.eqtls.MAF05.overlap.varsort.",chri,".txt"),
                col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
})

## write vars only to disk
overlap_vars_05 = data.frame(var = as.character(unique(cav_eqtls_05_var_sorted$var)), 
                          stringsAsFactors = FALSE)
write.table(overlap_vars_05, file="variant_sets/overlap_vars.MAF05.list",
            row.names=FALSE, quote=FALSE, col.names=FALSE)
## write vars + pos by chr
overlap_vars_05$chr = sapply(overlap_vars_05$var, function(var) {
  strsplit(var, '_')[[1]][1]
})
overlap_vars_05$pos = sapply(overlap_vars_05$var, function(var) {
  strsplit(var, '_')[[1]][2]
})
sapply(paste0('chr',seq(1,22)), function(chri) { 
  overlap_vars_05 %>% 
    filter(chr==chri) %>%
    select(chr, pos, var) %>%
    write.table(paste0("variant_sets/by_chr/overlap_vars.MAF05.",chri,".txt",
                       row.names=FALSE, quote=FALSE, col.names=TRUE))
})


nrow(cav_eqtls_05)
length(unique(cav_eqtls_05$var))
length(unique(cav_eqtls_05$gene))
tail(sort(table(cav_eqtls_05$gene)))
hist(table(cav_eqtls_05$gene))



cav_eqtls_10 = merge(caviar_eqtls, cav_overlap_10,
                     by.x=c('V2','V3'), by.y=c('V1','V2'))
names(cav_eqtls_10) <- c('chr','pos','gene','var')
cav_eqtls_10$chr <- factor(cav_eqtls_10$chr,
                           levels=paste0('chr',seq(1,22)))
cav_eqtls_10$pos <- as.integer(cav_eqtls_10$pos)
cav_eqtls_10_gene_sorted = cav_eqtls_10 %>%
  arrange(gene,chr,pos) %>%
  select(gene, chr, pos, var)
write.table(cav_eqtls_10_gene_sorted, "variant_sets/caviar_var.95set.eqtls.MAF10.overlap.genesort.txt",
            col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
sapply(paste0('chr',seq(1,22)), function(chri) { 
  cav_eqtls_10_gene_sorted %>% filter(chr==chri) %>%
    write.table(paste0("variant_sets/by_chr/caviar_var.95set.eqtls.MAF10.overlap.genesort.",chri,".txt"),
                col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
})
cav_eqtls_10_var_sorted = cav_eqtls_10 %>%
  arrange(chr,pos, gene) %>%
  select(chr, pos, gene, var)
write.table(cav_eqtls_10_var_sorted, "variant_sets/by_chr/caviar_var.95set.eqtls.MAF10.overlap.varsort.txt",
            col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
sapply(paste0('chr',seq(1,22)), function(chri) { 
  cav_eqtls_10_var_sorted %>% filter(chr==chri) %>%
    write.table(paste0("variant_sets/by_chr/caviar_var.95set.eqtls.MAF10.overlap.varsort.",chri,".txt"),
                col.names=TRUE, sep='\t', row.names=FALSE, quote=FALSE)
})

## write vars only to disk
overlap_vars_10 = data.frame(var = as.character(unique(cav_eqtls_10_var_sorted$var)), 
                          stringsAsFactors = FALSE)
write.table(overlap_vars_10, file="variant_sets/overlap_vars.MAF10.list",
            row.names=FALSE, quote=FALSE, col.names=FALSE)
## write vars + pos by chr
overlap_vars_10$chr = sapply(overlap_vars_10$var, function(var) {
  strsplit(var, '_')[[1]][1]
})
overlap_vars_10$pos = sapply(overlap_vars_10$var, function(var) {
  strsplit(var, '_')[[1]][2]
})
sapply(paste0('chr',seq(1,22)), function(chri) { 
  overlap_vars_10 %>% 
    filter(chr==chri) %>%
    select(chr, pos, var) %>%
    write.table(paste0("variant_sets/by_chr/overlap_vars.MAF10.",chri,".txt",
                     row.names=FALSE, quote=FALSE, col.names=TRUE))
})

nrow(cav_eqtls_10)
length(unique(cav_eqtls_10$var))
length(unique(cav_eqtls_10$gene))
tail(sort(table(cav_eqtls_10$gene)))
hist(table(cav_eqtls_10$gene))




