#!/bin/bash

#=============================================================
#
# Job script for finding single orthologues with OrthoFinder
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
#SBATCH --time=168:00:00
# Job name
#SBATCH --job-name=orthofinder-strep-test
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

orthofinder -f ../input -M msa -os -a 40 -t $SLURM_NTASKS

conda deactivate

#=========================================================
# Epilogue script to record job endtime and runtime
#=========================================================
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
