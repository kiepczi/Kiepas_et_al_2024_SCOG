import pandas as pd
from collections import Counter
from Bio import SeqIO
from pathlib import Path


datadir = Path("../output/SCOGs_sequences").expanduser()
filenames = sorted(datadir.glob("*"))

unique_nucleotide = {}

for _ in filenames:
    OG = _.stem
    records = list(SeqIO.parse(_, "fasta"))
    unique_nucleotide[OG] = len(set([_.seq for _ in records]))

count = Counter([v for k, v in unique_nucleotide.items()])

df = pd.DataFrame(list(unique_nucleotide.items()), columns=['orthogroup','unique_nucleotide_variants'])

print(df)


# datadir = Path("../output/sco_protein_trimmed_alignments").expanduser()
# filenames = sorted(datadir.glob("*"))


# unique_protein = {}

# for _ in filenames:
#     OG = _.stem.split('_')[0]
#     records = list(SeqIO.parse(_, "fasta"))
#     unique_protein[OG] = len(set([_.seq for _ in records]))


# df['unique_protein_variants'] = df['orthogroup'].map(unique_protein)

# print(unique_protein)

df.to_csv("SCOG_nucleotide_variant.csv", index=False)

# data = pd.read_csv("../output/SCOGs_distribution_vizualisation_data.csv")

# unique_nucleotide = {}

# for _ in data: 
#     if _ != "Name":
#         print(_, len(data[_]), len(set(data[_])))
#         unique_nucleotide[_] = len(set(data[_]))

# print(unique_nucleotide)

# print(Counter([v for k, v in unique_nucleotide.items()]))