"""This scrip was used to calculate the GC% content of all representative genomes, and the single copy orthogroups sequences.
"""


#Set Up
from Bio import SeqIO
import pandas as pd
from collections import Counter
from pathlib import Path


#Function that will calculate GC% content
def GC_content(sequence):
    """Return GC content for a given DNA sequence
    """

    GC_percentage = (sequence.count('C')+sequence.count('G'))/len(sequence)*100

    return GC_percentage


#Path to genomes gbk files. Here, we will you the input file for antismash analysis, as it has all gbk files for all represemtatve genomes used in this analysis
datadir = Path("../../supplementary_file_6/input")
filenames = sorted(datadir.glob("*"))

#Calculate genomes GC% content
genome_GC = {}

for _ in filenames:
    accession = '_'.join(str(_).split('/')[-1].split('_')[:2])
    genome_record = ''.join([str(_.seq) for _ in list(SeqIO.parse(_, 'gb'))])
    genome_GC[accession] = GC_content(genome_record)


df = pd.DataFrame(['_'.join(str(_).split('/')[-1].split('_')[:2]) for _ in filenames], 
                                              columns =['accession'])

#Assign this information to dataframe
df['genome_GC'] = df['accession'].map(genome_GC) 



datadir = Path("../../supplementary_file_5/output/sequences/nucleotide")
filenames = sorted(datadir.glob("*"))

for _ in filenames:
    if '.DS_Store' not in str(_):
        OG = str(_).split('/')[-1].split('_')[0]
        data = {_.description:GC_content(_.seq) for _ in list(SeqIO.parse(_, 'fasta'))}
        df[OG] = df['accession'].map(data)


#Melt data. Here, we melt data to create a tidy data
df = pd.melt(df, id_vars=["accession", 'genome_GC'], 
             value_name="SCO_GC").sort_values('SCO_GC')

#Renaming columns
df = df.rename(columns={'variable': 'Orthogroup'})

#Add information about the monophyly of the SCOGs variant 
monophyly = pd.read_csv("../../supplementary_file_10/output/SCOGs_distribution_vizualisation_data_reduced.csv")
monophyly = pd.melt(monophyly, id_vars=["Name"], 
             value_name="monophyletic").sort_values('monophyletic')
monophyly = monophyly.rename(columns={'variable': 'Orthogroup', 'Name':'accession'})

df = df.merge(monophyly, on=['accession', 'Orthogroup'], how='outer').fillna('Yes')



#Changing the values in the monophyletic_status column. If #F0F0F0 then change to 'Yes', else change to 'No'
def correct_monophyly(row):

    if row['monophyletic'] == '#F0F0F0':
        status = 'Yes'

    elif row['monophyletic'] == 'Yes':
        status = 'Yes'
    else:
        status = 'No'
    return status

# # Apply the custom function to 'Column1'
df['monophyletic'] = df.apply(correct_monophyly, axis=1)


df.to_csv('../output/GC_content.csv', index=False)