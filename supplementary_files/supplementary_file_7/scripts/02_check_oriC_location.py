# """This script was used to check the location of oriC for genomes assembled to complete or chromosomal level in NCBI.
# """

import pandas as pd
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


# #Loading data which contains information about the sequence type location of the SCOG. eg. here we are only interested in SCOGS present on hromosome, not scaffold or plasmid
# df = pd.read_csv("../output/SCOGs_general_info.csv")

# #Here, we only retain rows where the SCO is found on the chromosome
# df = df[df['seq_type'] == 'chromosome']

# #Now we can get a list of genomes of interest
# genomes = list(set([_ for _ in df['accession']]))

# #Now we can get oriC location for every genome of interest
# #Path to GenBank files
# datadir = Path("../../supplementary_file_1/NCBI_genomes")
# filenames = sorted(datadir.glob("strep*/GCF*/*.gbff"))


# oriC_location = {} #Hold dictionary to append the location of the oriC
# oriC_strand = {}
# for _ in filenames:
#     if str(_).split('/')[-2] in genomes:
#         records = list(SeqIO.parse(_, "genbank"))
#         chrom = records[0]
#         reverse_chrom = chrom.reverse_complement()
#         for feature in chrom.features:
#             if feature.type == "gene":
#                 if "gene" in feature.qualifiers and feature.qualifiers["gene"][0] == "dnaA":
#                     if feature.location.__dict__['_strand'] == 1:
#                         oriC_location[str(_).split('/')[-2]] = (feature.location.__dict__['_start']*100)/len(chrom)
#                         oriC_strand[str(_).split('/')[-2]] = 'positive'
#                     elif feature.location.__dict__['_strand'] == -1:
#                         oriC_location[str(_).split('/')[-2]] = (feature.location.__dict__['_start']*100)/len(chrom)
#                         oriC_strand[str(_).split('/')[-2]] = 'negative'

            
# df = pd.DataFrame()
# df['accession'] = genomes
# df['oriC_location'] = df['accession'].map(oriC_location) 
# df['oriC_strand'] = df['accession'].map(oriC_strand) 

# df.to_csv('../output/oriC_location.csv')

df = pd.read_csv("../output/oriC_location.csv")

figure(figsize=(15,20), dpi=80)

plt.scatter([_ for _ in df['oriC_location']], [_ for _ in df['accession']])
plt.xlabel('Chromosome length (%)', fontsize=14)
plt.ylabel('NCBI genome accession', fontsize=14)
# plt.title('oriC location')
plt.xlim(-1, 100)
# Increase the font size of the tick labels
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.savefig('../../supplementary_file_8.pdf', format='pdf')


