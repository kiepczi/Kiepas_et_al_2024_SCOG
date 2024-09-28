"""This script was used to determine the loction of SCOGs on the chromosome.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from collections import Counter

#Getting a list of genomes of interest
"""We have checked for the oriC location on the genomes. 3 out of 63 genomes investigated had oriC located on either arm of the chromosome, 
insread of in the centre. We did not find any publications stating that these were circular chromosomes, therefore were discareded from futher analysis.
"""

genomes_to_exclude = ["GCF_018128905.1", "GCF_009834125.1", "GCF_000147815.2"]
genomes = [genome for genome in pd.read_csv("../output/oriC_location.csv")['accession'] if genome not in genomes_to_exclude]


#Now we can get the SCOGs information. eg. the protein id
df = pd.read_csv("../output/SCOGs_general_info.csv")

#Here, we only retain rows where the SCO is found on the chromosome, and oriC is in the center (accession found in list of genomes)
df = df[df['accession'].isin(genomes)]

#To thid df we can now add the strand type on which the oriC is located. (eg. positive or negative)
df['oriC_strand'] = df['accession'].map(pd.read_csv('../output/oriC_location.csv').set_index('accession').to_dict()['oriC_strand'])

#Generating function that takes a GenBank file, checks if the gene 'dnaA' is on a reverse strand, and if so, reverse the genomic sequence, and save it as fna file.
datadir = Path("../../supplementary_file_1/NCBI_genomes/").expanduser()
filenames = sorted(datadir.glob("strep*/GCF*/*.gbff"))

def get_SCOG_location(row):
    """Return protein product."""

    file = ''
    for _ in filenames:
        if str(_).split('/')[-2] == row['accession']:
            file = str(_)

    records = list(SeqIO.parse(file, "genbank"))


    try:
        chrom = records[0]
        reverse_chrom = chrom.reverse_complement()
        if row['oriC_strand'] == 'positive':
            for feature in chrom.features:
                if feature.type == "CDS":
                    if "protein_id" in feature.qualifiers and feature.qualifiers["protein_id"][0] == row["protein_id"]:
                        SCO_loc = (feature.location.__dict__['_start']*100)/len(chrom)
                        print(len(chrom), row['protein_id'], file)
        elif row['oriC_strand'] =='negative':
            for feature in reverse_chrom.features:
                if feature.type == "CDS":
                    if "protein_id" in feature.qualifiers and feature.qualifiers["protein_id"][0] == row["protein_id"]:
                        SCO_loc = (feature.location.__dict__['_start']*100)/len(chrom)
                        print(len(chrom), row['protein_id'], file)
    except IndexError:
        SCO_loc = 'NA'


    return SCO_loc




df['SCOG_location'] = df.apply(get_SCOG_location, axis=1)



#Getting labels for visual representation
current_annotations = defaultdict(list)

for index, row in df.iterrows():
    current_annotations[row['Orthogroup']].append(row['product'])

labels = {}

for OG, annotations in current_annotations.items():
    labels[OG] = OG +': ' +str(Counter(annotations)).strip("Counter({'").replace("'})", "").replace("': ", '[').replace(", '", '] | ').replace('})', ']')

df['label'] = df['Orthogroup'].map(labels)

df.to_csv('../output/SCOGs_location.csv', index=False)




                        


