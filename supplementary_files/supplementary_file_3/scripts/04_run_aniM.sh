#!/bin/bash

#=============================================================
#
# Job script for running pyANI v3
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
#SBATCH --time=48:00:00
# Job name
#SBATCH --job-name=pyani_test
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

# Enable conda enviroment
module load anaconda/python-3.8.8/2021.05
source ~/.bashrc

#Orthofinder only works in orthofinder conda environment
conda activate pyani_v3
#=========================================================
# Prologue script to record job details
#=========================================================
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------

pyani anim -i ../input/ -o ../output/ -v -l ../output/pyani_representative_genomes.log --name pyani_295 --labels ../input/custom_labels.txt --classes ../input/custom_classes.txt --workers 40
conda deactivate

#=========================================================
# Epilogue script to record job endtime and runtime
#=========================================================
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
