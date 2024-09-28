#!/bin/bash

#=============================================================
#
# Job script for phylogenetic analysis
#
#=============================================================

#======================================================
# Propagate environment variables to the compute node
#SBATCH --export=ALL
# Project Account
#SBATCH --account=pritchard-grxiv
# Run in the standard partition (queue)
#SBATCH --partition=standard
# No of cores required (max. of 40)
#SBATCH --ntasks=40 --nodes=1
# Runtime (hard)
#SBATCH --time=60:00:00
# Job name
#SBATCH --job-name=MLSA
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

# Enable conda enviroment
module load anaconda/python-3.8.8/2021.05
source ~/.bashrc

#Orthofinder only works in orthofinder conda environment
conda activate orthofinder
#=========================================================
# Prologue script to record job details
#=========================================================
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------

raxml-ng --check --msa ../output/alignments/concatenated/no_gaps_concatenated_sco.fasta --model ../output/modeltest_ng_output/output.part.bic --prefix ../output/tree/01_check
raxml-ng --parse --msa ../output/alignments/concatenated/no_gaps_concatenated_sco.fasta --model ../output/modeltest_ng_output/output.part.bic --prefix ../output/tree/02_parse
raxml-ng --msa ../output/tree/02_parse.raxml.rba --model ../output/modeltest_ng_output/output.part.bic --all --bs-trees 100 --seed 1655486274 --prefix ../output/tree/03_bootstrap --threads $SLURM_NTASKS
raxml-ng --support --tree ../output/tree/03_bootstrap.raxml.bestTree --bs-trees ../output/tree/03_bootstrap.raxml.bootstraps --prefix ../output/tree/04_tbe --bs-metric tbe

conda deactivate

#=========================================================
# Epilogue script to record job endtime and runtime
#=========================================================
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
