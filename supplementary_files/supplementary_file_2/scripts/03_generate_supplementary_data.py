"""This script was used to generate supplementary data containing all information about the
genomes used in this manuscript.
"""

import pandas as pd
from collections import defaultdict
from collections import Counter
from pathlib import Path
from Bio import SeqIO


#Loading MLST data
MLST_data = pd.read_csv("../../supplementary_file_1/MLST_data/Genome_ST_info.csv")


#Loading 
df = pd.read_csv('../output/quality.txt', delim_whitespace=True, skiprows=3,  skipfooter=1,  header=None)
df.columns = [
        'Bin_ID', 'Marker', 'Lineage', 'genomes', 'markers', 'marker_sets', '0', '1', '2', '3', '4', '5+',
        'completeness', 'contamination', 'strain_heterogeneity'
    ]
#Get_accessions
accessions = {_:'_'.join(str(_).split('_')[0:2]) for _ in df['Bin_ID']}
df['accession'] = df['Bin_ID'].map(accessions)

MLST_data['completeness'] = MLST_data['accession'].map(df.set_index('accession')['completeness'].to_dict())
MLST_data['contamination'] = MLST_data['accession'].map(df.set_index('accession')['contamination'].to_dict())
MLST_data['strain_heterogeneity'] = MLST_data['accession'].map(df.set_index('accession')['strain_heterogeneity'].to_dict())

columns_to_delete = ['pyani_label', 'degree']

for _ in MLST_data:
    if _ in columns_to_delete:
        del MLST_data[_]

#Getting list of genomes of interest
genomes = [_ for _ in MLST_data['accession']]

#Path to genomes
datadir = Path("../../supplementary_file_1/NCBI_genomes/")
filenames = sorted(datadir.glob("streptom*/*/*feature_table.txt"))

total_genes = {}
total_cds = {}

for _ in filenames:
    genome = str(_).split('/')[-2]
    if genome in genomes:
        data = pd.read_csv(_, sep='\t')
        total_genes[genome] = Counter([_ for _ in data["# feature"] if _ =='gene'])['gene']
        total_cds[genome] = Counter([_ for _ in data["# feature"] if _ =='CDS'])['gene']

#Mapping values
MLST_data['total_genes'] = MLST_data['accession'].map(total_genes)    
MLST_data['total_cds'] = MLST_data['accession'].map(total_cds)    



# Calculating GC content
#Function that will calculate GC% content
def GC_content(sequence):
    """Return GC content for a given DNA sequence
    """

    GC_percentage = (sequence.count('C')+sequence.count('G'))/len(sequence)*100

    return GC_percentage


#Path to genomes
datadir = Path("../../supplementary_file_1/NCBI_genomes/")
filenames = sorted(datadir.glob("streptom*/*/*gbff"))

#Calculate genomes GC% content

genome_GC = {}
genome_size = {}

for file in filenames:
    accession = '_'.join(str(file).split('/')[-1].split('_')[:2])
    if accession in genomes:
        genome_data = ''.join([str(_.seq) for _ in list(SeqIO.parse(file, "genbank"))])
        genome_GC[accession] = GC_content(genome_data)
        genome_size[accession] = len(genome_data)

MLST_data['GC_content'] = MLST_data['accession'].map(genome_GC)    
MLST_data['genome_size'] = MLST_data['accession'].map(genome_size)    


MLST_data.to_csv("all_genomes_additional_data.csv", index=False)