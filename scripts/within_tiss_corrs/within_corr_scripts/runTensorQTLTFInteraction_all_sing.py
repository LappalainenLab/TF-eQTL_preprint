#!/usr/bin/env python3

## NOTES FOR SINGULARITY
## Singularity loaded using: -B 
##	~/.singularity/tensorqtl,
##	/gpfs/commons/home/eflynn/projects/TFi-eQTL:/tfi,
##	/gpfs/commons/home/eflynn/projects/TFi-eQTL_pilot2:/tfi_old,
##	/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/:/gtex,
## Script updated to include those prefixes



############################ Python 3.7.3, GTEx - INTERACTION

import os
import sys
import argparse

import numpy as np
import pandas as pd

import tensorqtl
from tensorqtl import genotypeio, cis, eigenmt

# def main():
parser = argparse.ArgumentParser()
parser.add_argument(
    '-t',
    '--tissue',
    dest='tissue',
    type=str,
    required=True,
    default=None,
    metavar='tissue',
    help="Tissue"
)
parser.add_argument(
    '-y',
    '--tf',
    dest='tf',
    type=str,
    required=True,
    default=None,
    metavar='tf',
    help="Transcription Factor (gene name description)"
)
if not sys.argv[1:]:
    sys.exit(parser.print_help())

args = vars(parser.parse_args()) # type: Dict[str, Any]
print(args)
print('%(tissue)s.%(tf)s.norm' % args)

prefix = '%(tissue)s.%(tf)s.norm.ieqtl.all_vars' % args
output_dir='/tfi/tensorqtl/%(tissue)s/%(tf)s' % args 

phenotype_bed_file='/tfi_old/phenotypes/normalized_expression/%(tissue)s.v8.normalized_expression.nochrX.bed.gz' % args
covariates_file='/gtex/eqtl/GTEx_Analysis_v8_eQTL_covariates/%(tissue)s.v8.covariates.txt' % args
interaction_file='/tfi_old/TF_expr/TF_norm_interaction_terms/%(tissue)s.norm.%(tf)s.txt' % args
genotype_plink_prefix='/tfi/genotypes/GTEx_v8.overlap_vars.MAF05'
#genotype_plink_prefix='/plink/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01'


## PHENOTYPES and COVARIATES
# Load phenotypes and covariates:
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T # samples x covariates
assert np.all(phenotype_df.columns==covariates_df.index)


## INTERACTION
# Load interaction data
interaction_s = pd.read_csv(interaction_file, sep='\t', index_col=0, header=0, squeeze=True)

## Select individuals that are in the interaction dataset
#phenotype_df = phenotype_df.iloc[:, phenotype_df.columns.isin(interaction_s.index)]
#covariates_df = covariates_df[covariates_df.index.isin(interaction_s.index)]
#assert np.all(phenotype_df.columns==covariates_df.index)

assert covariates_df.index.isin(interaction_s.index).all()
interaction_s = interaction_s.loc[covariates_df.index].astype(np.float32)


## GENOTYPES
# Genotypes can be loaded as follows, where plink_prefix_path is the path to the VCF in PLINK format:
# for VCFs with hard GT calls only, specify type as np.int8 to save memory:
## for all genotypes:
pr = genotypeio.PlinkReader(genotype_plink_prefix, select_samples=phenotype_df.columns) #, dtype=np.int8) 
## bug in the code - need to select the samples

## load into dataframes:
#genotype_df = pd.DataFrame(pr.get_all_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom','pos']]



## RUN TENSORQTL
## see code at https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/cis.py
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix,
	interaction_s=interaction_s, maf_threshold_interaction=0.05, 
	run_eigenmt=True, output_dir=output_dir) 

# if __name__ == '__main__':
#     main()
