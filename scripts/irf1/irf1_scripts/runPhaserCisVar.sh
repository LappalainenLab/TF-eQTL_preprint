#!/bin/bash
echo starting script
source ~/.bashrc
module unload python
module load python/2.7.8
cd ~/projects/ase_scripps/ase_rerun/
python phaser_cis_var.py --bed /gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/'/ase/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_matrix_WASP.gw_phased.txt.gz' --vcf /gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8//ase/GTEx_Analysis_v8_phASER/phASER_GTEx_v8_merged.vcf.gz --pairs test_pairs_topeQTLs.txt --map sample_maps/Colon_Sigmoid.txt --o results_phaser_cis_var/Colon_Sigmoid.txt --t 2 --ignore_v 1
echo Done in bash.
