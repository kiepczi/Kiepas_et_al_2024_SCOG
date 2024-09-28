"""This scipt was used to extract nucleotide sequences for all SCO needed for backtranslation. 
"""

#SetUp
from Bio import SeqIO
from pathlib import Path
from Bio import SeqIO
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pandas as pd
import os
from collections import Counter

#Getting a list of Single Copy Orthologues
#Load Single Copy Orthologue orthofinder output
with open('../../supplementary_file_4/output/Orthogroups/Orthogroups_SingleCopyOrthologues.txt') as f:
    SCOGs = [_.strip('\n') for _ in f.readlines()]


#Loading Orthogroup data (Orthofinder Output) and getting genome accessions as columns
SCO_data = pd.read_csv('../../supplementary_file_4/output/Orthogroups/Orthogroups.tsv', sep='\t')
new_data_cols = {_:'_'.join(_.split('_')[0:2]) for _ in SCO_data if _ != 'Orthogroup'}
SCO_data = SCO_data.rename(columns = new_data_cols)

#Slice data to only keep rows that represent SCOGs
SCO_data = SCO_data[SCO_data['Orthogroup'].isin(SCOGs)]

# Melt the DataFrame to have column names in a separate column
SCO_data = pd.melt(SCO_data, id_vars=['Orthogroup'], var_name='assembly_accession', value_name='protein_id')



#Extracting sequences
def extract_nucleotide_and_protein_seq(accession, protein_id):
    """Extract nucleotide and protein sequences.
    """

    print(accession)
    #Getting Protein Sequences 
    datadir = Path("../../supplementary_file_1/NCBI_genomes/streptomyces_faa_files")
    filenames_protein = [_ for _ in sorted(datadir.glob("*")) if accession in str(_)][0]
    prot_records = [SeqRecord(_.seq, description=accession, id=accession,name=accession) for _ in SeqIO.parse(filenames_protein, "fasta") if protein_id in _.description]


    #Getting Nucleotide Sequence and trimming stop codons
    stop_codons = ["TAG", "TAA", "TGA"]

    datadir = Path("../../supplementary_file_1/NCBI_genomes/")
    filenames_nucleotide = [_ for _ in sorted(datadir.glob("strep*/*/*cds_from_genomic.fna")) if accession in str(_)][0]
    nucleotide_records = [SeqRecord([_.seq[:-3] if _.seq[-3:] in stop_codons else _.seq][0], description=accession, id=accession,name=accession) for _ in SeqIO.parse(filenames_nucleotide, "fasta") if protein_id in _.description]


    return prot_records, nucleotide_records

SCOGs_protein_seq = defaultdict(list)
SCOGs_nucleotide_seq = defaultdict(list)

for index, row in SCO_data.iterrows():
    orthogroup = row['Orthogroup']
    accession = row['assembly_accession']
    protein_id = row['protein_id']
    protein_sequence, nucleotide_sequence = extract_nucleotide_and_protein_seq(accession, protein_id)
    SCOGs_protein_seq[orthogroup].extend(protein_sequence)
    SCOGs_nucleotide_seq[orthogroup].extend(nucleotide_sequence)



#Check if they are any duplicate sequneces, and if the sequences are the same. Save this to csv,
def are_all_items_same(input_list):
    # Check if the list is empty
    if not input_list:
        return "Yes"
    
    # Compare all items to the first item
    first_item = input_list[0]
    return "Yes" if all(item == first_item for item in input_list) else "No"

# Create an empty DataFrame
SCO_data_troublesome = {
    'Orthogroup': [],
    'assembly_accession': [],
    'copies': [],
    'protein_seq_identical': [],
    'nucleotide_seq_identical': []
}

SCO_data_troublesome = pd.DataFrame(SCO_data_troublesome)



for SCOG_group, nt_records in SCOGs_nucleotide_seq.items():
    for SCO_group, prot_records in SCOGs_nucleotide_seq.items():
        if SCOG_group == SCO_group:
                duplicates = list(set([accession for accession, count in Counter([_.description for _ in nt_records]).items() if count >=2]))
                for accession in duplicates:
                    troublesome_records_nt = [record for record in nt_records if record.description == accession]
                    troublesome_records_prot = [record for record in prot_records if record.description == accession]
                    if len(troublesome_records_nt) !=0:
                        SCO_data_troublesome.loc[len(SCO_data_troublesome)] = [SCOG_group, accession, len(nt_records), are_all_items_same([_.seq for _ in troublesome_records_prot]), are_all_items_same([_.seq for _ in troublesome_records_nt])]

# #Saving data with copies per each to CSV file
SCO_data_troublesome.to_csv("../output/additional_information/copies.csv", index=False)


#Getting a function that will check if the protein and nucleotide sequences have to be trimmed.
def trim_seq(prot, nt):
    """Trim protein and nucleotide sequences,
    so that the length of nucleotide sequence
    is exactly the 3x the length of the protein
    """
    if len(nt) - len(prot*3) < 0:
        prot_idx = len(nt) - len(prot*3)
        prot = prot[:prot_idx]
        nt_idx = len(prot*3) - len(nt)
        nt = nt[:nt_idx]
    elif len(nt) - len(prot*3) >= 1:
        nt_idx = len(prot*3)- len(nt)
        nt = nt[:nt_idx]
    else:
        nt = nt
        prot = prot

    return prot, nt



SCOGs_protein_seqs_checked = defaultdict(list)
SCOGs_nucleotide_seqs_checked = defaultdict(list)




for SCOG_group, nt_records in SCOGs_nucleotide_seq.items():
    for SCO_group, prot_records in SCOGs_protein_seq.items():
        if SCOG_group == SCO_group:
            seen = []
            for nt_record in nt_records:
                for prot_record in prot_records:
                    if nt_record.description == prot_record.description:
                        if nt_record.description not in seen:
                            fixed_prot, fixed_nt = trim_seq(prot_record.seq, nt_record.seq)
                            SCOGs_protein_seqs_checked[SCO_group].append(SeqRecord(fixed_prot, description=nt_record.description, id=nt_record.description,name=nt_record.description))
                            SCOGs_nucleotide_seqs_checked[SCO_group].append(SeqRecord(fixed_nt, description=nt_record.description, id=nt_record.description,name=nt_record.description))
                            seen.append(nt_record.description)

for OG_group, records in SCOGs_nucleotide_seqs_checked.items():
     SeqIO.write(records, f'../output/sequences/nucleotide/{OG_group}_nt_trimmed.fasta', "fasta")

for OG_group, records in SCOGs_protein_seqs_checked.items():
     SeqIO.write(records, f'../output/sequences/protein/{OG_group}_prot_trimmed.fasta', "fasta")

