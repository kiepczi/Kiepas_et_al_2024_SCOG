"""This script was used to concatenate CDS sequences, and generate a single 
aligmnents file that can be used as an output for RAXML-ng. Additionally, this
scipt was used to generate partitional file for Modeltest-ng. 

The partition file format is as follow:

DATA_TYPE, PARTITION_NAME = PARTITION_SITES
"""

#SetUp
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Path to SCO cds alignments
datadir = Path("../output/alignments/SCO_backthreaded")
filenames = sorted(datadir.glob("OG*"))


#Concatenate sequence and generate partition file
concatenated = {}
partitions = []
last_position = 0
for filename in filenames:
    OG = str(filename).split('/')[-1].split(' ')[0]
    records = SeqIO.parse(filename, "fasta")
    for record in records:
        accession = record.description
        sequence = record.seq

        if accession not in concatenated:
            concatenated[accession] = sequence
        else:
            concatenated[accession] += sequence

    next_position = last_position + len(record)
    partitions.append(f"DNA, {OG} = {last_position +1}-{next_position}")
    last_position = next_position

#Save concatenate sequences
sequence_list = []
for key, value in concatenated.items():
    record = SeqRecord(
                seq = value, 
                description = key,
                id = key

            )
    sequence_list.append(record)
SeqIO.write(sequence_list, "../output/alignments/concatenated/initial_concatenated_sco.fasta", "fasta")

#Save partition file

with open("../output/alignments/concatenated/initial_concatenated_modeltest.part", "w") as text_file:
    text_file.write('\n'.join(partitions))

