#!/bin/bash
##
##	calcMeffGeneGTs_array_wrapper.sh
##
##	EDF 5/21/20
##

which Rscript

Rscript scripts/calcMeffGeneGTs_array.R $SLURM_ARRAY_TASK_ID

echo "Done."
date


