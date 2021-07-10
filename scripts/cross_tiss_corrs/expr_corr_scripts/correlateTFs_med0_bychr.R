#!/usr/bin/Rscript
##
##  correlateTFs_med0_bychr.R
##
##  Correlate TF expr with eqtl effect size.
##
##  USAGE: ./correlateTFs_med0_bychr.R <maf> <tf name> <chr number>
##                                  e.g. 05 ATF3 20
##
##  EDF 4/28/20
##

## Load libraries
library(dplyr)
library(seqminer)

## Set working directory
setwd("/gpfs/commons/groups/lappalainen_lab/eflynn/projects/crosstiss_tf_corrs/")

## Get command line arguments
args = commandArgs(trailingOnly=TRUE)
maf=args[1]
#maf="05"
tf=args[2]
#tf="FEZF1"
chr=args[3]
#chr=4

print(paste("maf is",maf))
print(paste("tf is",tf))
print(paste("chr is",chr))
print("med0 eGenes only")

print("Loading data...")

## Load tissue names table and tissue info table
tiss_trans = read.table("input_files/tissue_translation_colors_v8.txt",
                        header=TRUE,sep='\t') %>%
  merge( read.table("input_files/tissue_sample_count.txt",
                    header=FALSE, sep='\t'),
         by.x=c('TISSUE_NAME'), by.y=c('V1') ) %>%
  rename(num_samps = V2)

## Load eqtl effect sizes
aFCs = read.table(paste0("afcs/combined/by_chr/all_tiss.afcs.cov.overlap_vars.MAF",maf,".chr",chr,".txt"),
                  header=TRUE, sep='\t')

## Load all gene expression
expr = read.table("input_files/genes.rnaseqc.median_tpm.all_tissues_v8.txt.gz",
                  header=TRUE, sep='\t')

## Load TF info
tfs = read.table("input_files/TFs.info.curated.txt",
                 header=TRUE, sep='\t')
tf_chipseq_file = paste0("input_files/encode_overlap/overlaps/ENCODE_TF_ChIPseq_overlap.GTEx_v8.curated_set.MAF",maf,".txt.gz")
tf_motif_file = paste0("input_files/motif_overlap/overlaps/HOCOMOCO_TF_motif_overlap.either.GTEx_v8.curated_set.MAF",maf,".txt.gz")


## Tissues found in effect size file (fewer than in expr file)
tissues_to_select = names(aFCs)[5:53]

## TF gene name
tf_gene = tfs %>%
  filter( TF == tf ) %>%
  pull(Name) %>%
  as.character()

## TF expr values for relevant tissues
tf_expr = expr %>% 
  filter(gene == tf_gene) %>% 
  select(all_of(tissues_to_select)) %>%
  unlist()


print("Beginning correlations...")
chr_old="chr0"
i=1

## Cycle through aFC table rows and calculate correlation for each
## Data from each row (dataframe at end of apply statement) 
##   will be 'rbind' together to form a large table
tf_corrs = do.call( 'rbind', apply(aFCs, 1, function(aFC_row) {
  
  print(aFC_row[1:5])
  
  ## Get variant of eqtl row
  var_i=as.character(aFC_row[1])
  #print(var_i)
  
  ## Check if chr is changing (not really needed since this is by chr)
  chr_i=strsplit(var_i,"_")[[1]][1]
  if (chr_i != chr_old) {print(paste("  Now on",chr_i))}
  chr_old <<- chr_i
  
  ## Print row number if we are on 1000/2000/etc. row
  if (i %% 1000 == 0) {print(paste("    Now on line",i))}
  i <<- i+1
  
  ## Get more eqtl info from variant/eqtl row
  pos_i=strsplit(var_i,"_")[[1]][2]
  gene_i=as.character(aFC_row[4])
  gene_name_i = expr %>%
    filter(gene == gene_i) %>%
    pull(description) %>%
    as.character()
  #print(gene_i)
  
  ## Select aFCs only (not eqtl info)
  afcs_i=as.numeric(aFC_row[5:53])
  #print(afcs_i)
  
  ## Select eGene expression for relevant tissues
  exprs_i = expr %>% 
    filter(gene == gene_i) %>% 
    select(all_of(tissues_to_select)) %>%
    unlist()
  #print(exprs_i)
  
  ## Get number of tissues with values for both afc and TF expr
  num_tiss = sum(!is.na(afcs_i) & !is.na(tf_expr) & exprs_i > 0)
  
  ## Select tissues that have values for both afc and TF expr
  cols_to_keep = which(!is.na(afcs_i) & !is.na(tf_expr) & exprs_i > 0)
  tissues_to_keep = tissues_to_select[cols_to_keep]
  
  ## Select afcs, eGene expr, tf expr only for tissues above
  afcs_i = afcs_i[cols_to_keep]
  exprs_i = exprs_i[cols_to_keep]
  tf_exprs_i = tf_expr[tissues_to_keep]
  
  if (num_tiss > 0) {
    ## Calculate afc stats
    max_afc = afcs_i[which.max(abs(afcs_i))]
    med_afc = median(afcs_i)
    min_afc = afcs_i[which.min(afcs_i*sign(max_afc))]
    afc_range = abs(max_afc - min_afc)
    
    ## Calculate eGene expr stats
    min_expr = min(exprs_i)
    med_expr = median(exprs_i)
    max_expr = max(exprs_i)
    expr_lfc = log2((max_expr+0.001)/(min_expr+0.001))
    
    ## Calculate TF expr stats
    min_tf_expr = min(tf_exprs_i)
    med_tf_expr = median(tf_exprs_i)
    max_tf_expr = max(tf_exprs_i)
    tf_expr_lfc = log2((max_tf_expr+0.001)/(min_tf_expr+0.001))
  } else {
    max_afc = med_afc = min_afc = afc_range = NA
    min_expr = med_expr = max_expr = expr_lfc = NA
    min_tf_expr = med_tf_expr = max_tf_expr = tf_expr_lfc = NA
  }
  
  ## If there are at least 3 tissues
  if (num_tiss > 2) {
    
    ## Run spearman correlation on afc, TF expr
    sp_cor = cor.test(afcs_i, tf_exprs_i)
    
    ## Get direction of correlation based on afc values
    ##   (are afcs getting closer or farther from 0?)
    if (med_afc > 0 & min_afc > -0.5*max_afc) {
      cor_dir_temp = sp_cor$estimate
    } else if (med_afc < 0 & min_afc < -0.5*max_afc) {
      cor_dir_temp = sp_cor$estimate * -1
    } else { cor_dir_temp = 0 }
    cor_dir = if (is.na(cor_dir_temp)) {NA
      } else if (cor_dir_temp > 0) {'pos'
      } else if (cor_dir_temp < 0) {'neg'
      } else {'unc'}
    
  } else { 
    ## Set null values if we don't do the correlation
    sp_cor = data.frame(estimate = NA, p.value=NA); cor_dir = NA 
  }
  
  ## Get chipseq overlap for this variant/TF
  chip_var = tabix.read.table( tf_chipseq_file, paste0(chr_i,":",pos_i,"-",pos_i) ) %>%
    filter(VARIANT==var_i) 
  chip_tf = chip_var %>% pull(!!tf)
  chip_sum = chip_var %>% pull(sum)
  
  ## Get motif overlap for this variant/TF
  motif_var = tabix.read.table( tf_motif_file, paste0(chr_i,":",pos_i,"-",pos_i) ) %>%
    filter(snp==var_i) 
  motif_tf = motif_var %>% pull(!!tf)
  motif_sum = motif_var %>% pull(sum)
  
  ## Return a single row dataframe of our data
  data.frame(
    gene = gene_i, gene_name = gene_name_i, variant = var_i, 
    min_afc = min_afc, med_afc = med_afc, max_afc = max_afc, afc_range = abs(max_afc - min_afc),
    min_expr = min_expr, med_expr = med_expr, max_expr = max_expr, expr_lfc = expr_lfc,
    tf = tf, tf_gene = tf_gene,
    min_tf_expr = min_tf_expr, med_tf_expr = med_tf_expr, max_tf_expr = max_tf_expr, tf_expr_lfc = tf_expr_lfc,
    num_tiss = num_tiss,
    sp_rho = sp_cor$estimate, 
    sp_p = sp_cor$p.value,
    sp_dir = cor_dir,
    tf_chip = chip_tf, all_chip = chip_sum,
    tf_motif = motif_tf, all_motif = motif_sum
  )
}) )


## Write table to file
write.table(tf_corrs,
            paste0("correlations/by_tf/by_chr/cross_tiss_tf_expr_corrs.med0.curated_set.MAF",maf,".",tf,".chr",chr,".txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')


# 
# 
# ### quick plot. can use for debugging individual cases
# tf_expr_t = expr %>%
#   filter(gene == tf_gene) %>%
#   select(tissues_to_select)
# var_t = "chr20_2510104_T_C_b38"
# gene_t = "ENSG00000088876.11"
# gene_name_t = expr %>%
#   filter(gene == gene_t) %>%
#   pull(description) %>%
#   as.character()
# afcs_t=aFCs %>%
#   filter(sid==var_t, pid==gene_t) %>%
#   select(tissues_to_select)
# 
# cbind(t(tf_expr_t), t(afcs_t)) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="tiss") %>%
#   ggplot(aes(log10(V1), V2, col=tiss)) +
#   geom_hline(yintercept=0) +
#   geom_point() +
#   scale_color_manual(values=as.character(tiss_trans[order(tiss_trans$TISSUE_ABBRV),'TISSUE_RCOL'])) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   xlab(paste0("log10(",tf," tpm)")) +
#   ylab("aFC") +
#   ggtitle(paste(gene_t, gene_name_t, "\n", var_t))
# 
# exprs_t = expr %>%
#   filter(gene==gene_t) %>%
#   select(tissues_to_select)
# cbind(t(exprs_t), t(afcs_t)) %>%
#   as.data.frame() %>%
#   rownames_to_column(var="tiss") %>%
#   ggplot(aes(log10(V1), V2, col=tiss)) +
#   geom_hline(yintercept=0) +
#   geom_point() +
#   scale_color_manual(values=as.character(tiss_trans[order(tiss_trans$TISSUE_ABBRV),'TISSUE_RCOL'])) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   xlab(paste0("log10(",gene_name," tpm)")) +
#   ylab("aFC") +
#   ggtitle(paste(gene_t, gene_name_t,"\n", var_t))
