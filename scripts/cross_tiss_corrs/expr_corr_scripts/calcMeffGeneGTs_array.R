#!/usr/bin/Rscript
##
##  calcMeffGeneGTs.R
##
##  EDF 5/11/20
##

library(poolr)
library(dplyr)
library(seqminer)

setwd("~/projects/crosstiss_tf_corrs/")

count = as.numeric(commandArgs(trailingOnly = TRUE)[1])
print(paste("Array number is",count))


genotype_file = "afcs/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf.gz"

var_genes = read.table("variant_sets/caviar_var.95set.eqtls.MAF05.overlap.genesort.txt",
                       header=TRUE, sep='\t')
genes = as.character(levels(var_genes$gene))
genes = genes[ (count*322 - 321):(count*322) ]
genes = genes[!is.na(genes)]

gene_tests = do.call('rbind', lapply(genes, function(gene_i) {
  print(gene_i)
  var_gene_i = var_genes %>%
    filter(gene==gene_i)
  #print(var_gene_i)
  
  gts_i = do.call('rbind', apply(var_gene_i, 1, function(row) {
    #print(row)
    tabix.read.table(genotype_file, paste0(row[2],":",row[3],"-",row[3])) %>%
      filter(ID == row[4]) %>%
      select(-c(1,2,3,4,5,6,7,8,9)) %>%
      mutate_all(function(x) { sum(as.numeric(strsplit(x,"/")[[1]])) })
  } ) )
  
  n_var = nrow(gts_i)
  cor_mat = cor(t(gts_i), use="pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  ## Is this legitimate? Change all NA values to 0.
  ## Other methods use shrinkage approximations of matrices (?)
  
  meff = meff(cor_mat, method='gao', C=0.99)
  
  
  ## cor(t(gts_i), : generate correlation matrix of vars vs. vars
  ##    use="pairwise.complete.obs"), : for each pair of corrs, 
  ##          use all entries that are not NA in either
  ##          (instead of deleting everything in the whole matrix if one of the entries is NA)
  ## method='gao', : use the Gao method to compute meff
  ## C=0.99 : eigenvectors must explain at least 99% of variance
  
  data.frame(gene=gene_i, n_var=n_var, meff=meff)
  
} ) )

write.table(gene_tests, paste0("correlations/meff_gts_split/gene_tests.meff_gts_gao_99.",count,".txt"),
            col.names=FALSE, row.names=FALSE, 
            quote=FALSE, sep='\t')


print("Rscript done.")

