"""This script was used to geenrate vizualisation data for SCOG location 
on chromosome where some nucletide varaiants have been present multiple times.
"""

import pandas as pd

#Data with location of SCOGs on chromosome
location_data = pd.read_csv("../output/SCOGs_location.csv")

#Data with multiple variants but found to be monophyletic
monophyletic_SCOGs = pd.read_csv("../../supplementary_file_10/output/SCOGs_distribution_vizualisation_data_monophyletic.csv")
monophyletic_SCOGs = pd.melt(monophyletic_SCOGs, id_vars=['Name'], var_name='Orthogroup', value_name='tag')

monophyletic_SCOGs = monophyletic_SCOGs.loc[monophyletic_SCOGs['tag'] != '#F0F0F0']
monophyletic_SCOGs = monophyletic_SCOGs.rename(columns={'Name': 'accession'})
monophyletic_SCOGs["HGT_status"] = "monophyletic"
del monophyletic_SCOGs["tag"]

#Data with multiple variants but found to be non monophyletic, but within the same genus
non_monophyletic_SCOGs_same_genus = pd.read_csv("../../supplementary_file_10/output/SCOGs_distribution_vizualisation_data_same_genus.csv")
non_monophyletic_SCOGs_same_genus = pd.melt(non_monophyletic_SCOGs_same_genus, id_vars=['Name'], var_name='Orthogroup', value_name='tag')

non_monophyletic_SCOGs_same_genus = non_monophyletic_SCOGs_same_genus.loc[non_monophyletic_SCOGs_same_genus['tag'] != '#F0F0F0']
non_monophyletic_SCOGs_same_genus = non_monophyletic_SCOGs_same_genus.rename(columns={'Name': 'accession'})
non_monophyletic_SCOGs_same_genus["HGT_status"] = "non_monophyletic_same_genus"
del non_monophyletic_SCOGs_same_genus["tag"]


#Data with multiple variants but found to be non monophyletic, but within diffrent genus
non_monophyletic_SCOGs_diffrent_genus = pd.read_csv("../../supplementary_file_10/output/SCOGs_distribution_vizualisation_data_genus_split.csv")
non_monophyletic_SCOGs_diffrent_genus = pd.melt(non_monophyletic_SCOGs_diffrent_genus, id_vars=['Name'], var_name='Orthogroup', value_name='tag')

non_monophyletic_SCOGs_diffrent_genus = non_monophyletic_SCOGs_diffrent_genus.loc[non_monophyletic_SCOGs_diffrent_genus['tag'] != '#F0F0F0']
non_monophyletic_SCOGs_diffrent_genus = non_monophyletic_SCOGs_diffrent_genus.rename(columns={'Name': 'accession'})
non_monophyletic_SCOGs_diffrent_genus["HGT_status"] = "non_monophyletic_diffrent_genus"
del non_monophyletic_SCOGs_diffrent_genus["tag"]

#Combining all 
df_combined = pd.concat([non_monophyletic_SCOGs_diffrent_genus, non_monophyletic_SCOGs_same_genus, monophyletic_SCOGs], axis=0, ignore_index=True)

#Making location data
new_loc_data = pd.merge(df_combined, location_data, on=['Orthogroup', 'accession'], how='left')
new_loc_data = new_loc_data.dropna(subset=['SCOG_location'])

new_loc_data.to_csv("../output/SCOGs_location_HGT.csv", index=False)