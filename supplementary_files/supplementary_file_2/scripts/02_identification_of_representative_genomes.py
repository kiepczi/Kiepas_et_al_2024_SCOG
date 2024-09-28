"""This script was used to map genome quality information to MLST generated data. 
"""

import pandas as pd
from collections import defaultdict
from collections import Counter


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

check_completeness = defaultdict(list)

for index, row in MLST_data.iterrows():
    species = row['pyani_species_ID']
    completeness = row['contamination']
    check_completeness[species].append(completeness)



# Sort the DataFrame by 'completeness' (highest to lowest), 'contamination' (lowest to highest), and 'strain_heterogeneity' (lowest to highest)
MLST_data = MLST_data.sort_values(by=['completeness', 'contamination'], ascending=[False, True])


# Drop duplicates based on 'pyani_species_ID' while keeping the first occurrence (highest completeness, lowest contamination)
MLST_data = MLST_data.drop_duplicates(subset='pyani_species_ID', keep='first')

# Reset the index of the resulting DataFrame
MLST_data = MLST_data.reset_index(drop=True)

MLST_data.to_csv('../output/representative_genomes.csv', index=False)

