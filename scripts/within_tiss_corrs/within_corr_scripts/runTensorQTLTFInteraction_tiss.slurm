#!/bin/bash
#SBATCH --partition=gpu          				# Specify GPU partition (queues)
#SBATCH --mem=10gb              	                        # Job memory request (30 works for most)
#SBATCH --time=2:00:00						# Time limit 4 hours
## #SBATCH --output=log/tensorqtl_TFnormi_all_%A_%a.output 	# By default stderr and stdout is redirected to the same file, use %j if not using task array
#SBATCH --gres=gpu:tesla:1					# To request for 1 GPU card (Tesla V100 only for pe2dg2 nodes)
## #SBATCH --nodelist=pe2dg2-02					# Special condition so I don't hog the cluster
## #SBATCH --job-name=TFnormiTensorQTL_tiss			# Job name (specified in command line now)
## #SBATCH --array=1-3 ##446					# Array range (number of TFs) (specified in command line now)


## Run: sbatch runTensorQTLTFNormInteraction.slurm
## NB! Need to ssh onto the login node (pe2-login)


module load singularity
#singularity run --nv \
#	-B ~/.singularity/tensorqtl,/gpfs/commons/home/eflynn/projects/TFi-eQTL:/tfi,/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/:/gtex,/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/plink_genotypes/:/plink \
#	~/.singularity/tensorqtl/tensorqtl.sif


## Prepare folders
cd ~/projects/TFi-eQTL
tiss=$1
declare -a TFs=(`awk '{print $1}' input_files/TFs.info.curated.txt | tail -n +2`)
tf=${TFs[$(($SLURM_ARRAY_TASK_ID - 1))]}

echo Tissue is $tiss
echo Transcription factor is $tf
echo tensorqtl version is `singularity exec --nv ~/.singularity/tensorqtl/tensorqtl.sif head -n 1 /opt/tensorqtl/tensorqtl/__init__.py`

## run tiss tf ieqtls
if [[ ! -e tensorqtl/$tiss/$tf/$tiss.$tf.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz ]]
then
	if [[ ! -d tensorqtl/$tiss/$tf ]]
	then
		mkdir tensorqtl/$tiss/$tf
	elif [[ -e tensorqtl/$tiss/$tf/$tiss.$tf.norm.ieqtl.all_vars.cis_qtl_pairs.chr1.parquet ]]
	then
		echo removing old parquet files
		srm tensorqtl/$tiss/$tf/$tiss.$tf.norm.ieqtl.all_vars.cis_qtl_pairs.*parquet
	fi
	echo running $tiss interactions
	singularity exec --nv \
		-B ~/.singularity/tensorqtl,/gpfs/commons/home/eflynn/projects/TFi-eQTL:/tfi,/gpfs/commons/home/eflynn/projects/TFi-eQTL_pilot2/:/tfi_old,/gpfs/commons/datasets/controlled/GTEx/dbgap_restricted/data/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/:/gtex \
		~/.singularity/tensorqtl/tensorqtl.sif \
		python3 scripts/runTensorQTLTFInteraction_all_sing.py --tissue $tiss --tf $tf
else
	echo $tiss output already exists.
	ls -l tensorqtl/$tiss/$tf/$tiss.$tf.norm.ieqtl.all_vars.cis_qtl_top_assoc.txt.gz
fi


echo Done.
date
