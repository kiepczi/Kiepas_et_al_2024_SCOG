"""This scripts was used to get input genomes.
"""
#Set Up
import pandas as pd
from pathlib import Path
import shutil

#Getting list of genomes of interest
representative_genome_data = [_ for _ in pd.read_csv("../../supplementary_file_2/output/representative_genomes.csv")['accession']]


datadir = Path("../../supplementary_file_1/NCBI_genomes/streptomyces_faa_files")
filenames = sorted(datadir.glob("*faa"))

for fname in filenames:
    if '_'.join(str(fname).split('/')[-1].split('_')[0:2]) in representative_genome_data and 'from' not in str(fname):
        print(fname)
        shutil.copy(fname, "../input/")
