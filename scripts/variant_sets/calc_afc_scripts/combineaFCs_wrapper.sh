#!/bin/bash
##
##	combineaFCs_wrapper.sh
##
##	EDF 4/22/20
##

maf=$1
chr=$SLURM_ARRAY_TASK_ID

/nfs/sw/R/R-3.6.0/bin/Rscript scripts/combineaFCs.R $maf $chr

echo "Done."
date


