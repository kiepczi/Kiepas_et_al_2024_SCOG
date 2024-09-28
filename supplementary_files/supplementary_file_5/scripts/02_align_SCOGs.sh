#!/bin/bash
# align_sco_mafft.sh
# Align single-copy orthologue sequences using MAFFT

# Create output directory
mkdir -p ../output/alignments/sco_protein_trimmed_alignments

# Align each set of SCOs
for fname in ../output/sequences/protein/*; do
    mafft $fname > ../output/alignments/sco_protein_trimmed_alignments/`basename ${fname%%_prot_trimmed.fasta}`_aligned.fasta
done

