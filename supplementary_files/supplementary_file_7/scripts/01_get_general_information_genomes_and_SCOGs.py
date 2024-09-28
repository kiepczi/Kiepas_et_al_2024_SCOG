"""This script was used to generate data with information about each SCOG per genome. 
Here we will check:
- Protein IDs assigned in NCBI/gbk files
- Product Names
- Record IDs
- Type of sequence where the SCOG is located eg. chromosome, plasmid, contig or scaffold
"""


import pandas as pd
from pathlib import Path
from Bio import SeqIO

#Load orthofinder ouput
df = pd.read_csv('../../supplementary_file_4/output/Orthogroups/Orthogroups.tsv', sep='\t')

#Load Single Copy Orthologue orthofinder output
with open('../../supplementary_file_4/output/Orthogroups/Orthogroups_SingleCopyOrthologues.txt') as f:
    lines = f.readlines()

#Get list of orthogroups that are SCO
sco = [_.replace('\n', '') for _ in lines]
#Slice data to only keep rows that represent SCO
SCO_data = df[df['Orthogroup'].isin(sco)]


#Rename columns to only have the genome accession 
new_data_cols = {_:'GCF_' + _.split('_')[1] for _ in SCO_data if _ != 'Orthogroup'}
SCO_data = SCO_data.rename(columns = new_data_cols)


# Assuming your DataFrame is named df
SCO_data = pd.melt(SCO_data, id_vars=['Orthogroup'], var_name='accession', value_name='protein_id')


#Path to all CDS files
datadir = Path("../../supplementary_file_1/NCBI_genomes/")
filenames = sorted(datadir.glob("strep*/GCF*/GCF*fna"))

def get_protein_product(row, col):
    """Return protein product."""

    file = ''
    for _ in filenames:
        if str(_).split('/')[-2] == row['accession'] and 'cds_from' in str(_):
            file = str(_)

    fasta = list(SeqIO.parse(file, 'fasta'))

    print(file)


    try:
        temporary = [_.description for _ in fasta if row['protein_id'] in _.description]
        product_info = temporary[0].split('protein=')[-1].split(']')[0]
        record_id = temporary[0].split('lcl|')[-1].split('_cds')[0]
        locus_tag = temporary[0].split('locus_tag=')[-1].split(']')[0]

    except IndexError:
        product_info = 'NA'
        record_id = 'NA'
        locus_tag = 'NA'


    return product_info, record_id, locus_tag





SCO_data['product'], SCO_data['record_id'], SCO_data['locus_tag']  = zip(*SCO_data.apply(get_protein_product, col='protein_id', axis=1))





#Path to all NCBI feature tables
datadir = Path("../../supplementary_file_1/NCBI_genomes")
filenames = sorted(datadir.glob("strep*/GCF*/GCF*table.txt"))


def get_seq_type(row):
    """Return protein product."""

    file = ''
    for _ in filenames:
        if str(_).split('/')[-2] == row['accession']:
            file = str(_)

    df = pd.read_csv(file, sep='\t')


    try:
        seq_type = df[(df['product_accession'] == row['protein_id'])& (df['locus_tag'] == row['locus_tag']) & (df['genomic_accession'] == row['record_id'])]['seq_type'].values[0] if not df.empty else None

    except IndexError:
        seq_type = 'NA'


    return seq_type


SCO_data['seq_type'] = SCO_data.apply(get_seq_type, axis=1)

SCO_data.to_csv('../output/SCOGs_general_info.csv', index=False)

