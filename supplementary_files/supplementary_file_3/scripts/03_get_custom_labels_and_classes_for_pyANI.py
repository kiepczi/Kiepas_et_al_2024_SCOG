"""This script was used to customise classes and label input files for pyANI analysis.

Here, we specify the pyANI grpup ID as class, and the labels will consists of information like
organism name, strain and the assigned STs.
"""

#Set Up
import pandas as pd
from pathlib import Path

#Load data

df = pd.read_csv(Path("../../supplementary_file_2/output/representative_genomes.csv"))

#Create a dictionary with genome accession/id as keys, and pyani custom label as values. 
# This will allow to match the new labels, with the ones provided by `pyani index` parameter
pyani_labels = df.set_index('accession').to_dict()['pyani_label']
pyani_classes = df.set_index('accession').to_dict()['pyani_group_ID']

#Create custom_labels.txt file for pyani analysis
def custom_pyani_labels(file, new_pyani_labels):
    """Write a tab separated file for a custom/alternative pyani's
    graphical output.
    """
    #Load data
    data = pd.read_csv(file, sep='\t', names=['MD5_hash', 'genome_file', 'label'])
    #Get a dictionary with genome file name as key, and current pyani label
    pyani_labels = data.set_index('genome_file').to_dict()['label']

    #Change pyani_labels values (current pyani labels), to a custom label if genome accession (new_pyani_label - key) is present in genome file (pyani_labels key)
    for genome_file, current_label in pyani_labels.items():
        for genome_accession, new_label in new_pyani_labels.items():
            if genome_accession in genome_file:
                pyani_labels[genome_file] = new_label
    
    #change 'label' column to contain a new custom label by matching the key of pyani_labels dictionary with data['genome_file] column values
    data['label'] = data['genome_file'].map(pyani_labels)

    #Saving the new file 
    file_name = str(file).replace('labels.txt', 'custom_labels.txt').replace('classes.txt', 'custom_classes.txt')


    data.to_csv(file_name, sep='\t', header=False, index=False)

    return print(data)

# custom_pyani_labels(('pyani_other/venezuelae/classes.txt'))
#Get path to existing pyani labels
datadir = Path("../input/")
filenames = sorted(datadir.glob("labels.txt"))

for _ in filenames:
    custom_pyani_labels(_, pyani_labels)



datadir = Path("../input/")
filenames = sorted(datadir.glob("classes.txt"))

for _ in filenames:
    custom_pyani_labels(_, pyani_classes)