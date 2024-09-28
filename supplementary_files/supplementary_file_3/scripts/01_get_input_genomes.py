"""This scripts was used to get input genomes.
"""
#Set Up
import pandas as pd
from pathlib import Path
import shutil

#Getting list of genomes of interest
representative_genome_data = [_ for _ in pd.read_csv("../../supplementary_file_2/output/representative_genomes.csv")['accession']]

print(len(representative_genome_data))

datadir = Path("../../../Kiepas_et_al_2023_MLST/supplementary_files/supplementary_file_2/data")
filenames = sorted(datadir.glob("strep*/GCF*/GCF*genomic.fna"))

for fname in filenames:
    if str(fname).split('/')[-2] in representative_genome_data and 'from' not in str(fname):
        shutil.copy(fname, "../input/")
