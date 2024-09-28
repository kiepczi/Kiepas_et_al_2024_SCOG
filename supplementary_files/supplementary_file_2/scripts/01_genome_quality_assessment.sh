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
#SBATCH --time=167:00:00
# Job name
#SBATCH --job-name=orthofinder-strep-test
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

# Enable conda enviroment
module load anaconda/python-3.8.8/2021.05
source ~/.bashrc

#Orthofinder only works in orthofinder conda environment
conda activate myenv
#=========================================================
# Prologue script to record job details
#=========================================================
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------

checkm tree ../input ../output -t $SLURM_NTASKS
checkm tree_qa ../output -f ../output/marker_file.txt
checkm lineage_set ../output ../output/marker_file.txt
checkm analyze ../output/marker_file.txt ../input ../output -t $SLURM_NTASKS
checkm qa ../output/marker_file.txt ../output -f ../output/quality.txt -t $SLURM_NTASKS


conda deactivate

#=========================================================
# Epilogue script to record job endtime and runtime
#=========================================================
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
