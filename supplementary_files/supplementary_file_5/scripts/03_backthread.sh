#!/bin/bash
# backtranslate_alignments.sh
# Align single-copy orthologue sequences using MAFFT

# Create output directory
mkdir -p ../output/alignments/SCO_backthreaded

# Align each set of SCOs
for fname in ../output/sequences/nucleotide/*; do
    t_coffee -other_pg seq_reformat -in ${fname} -in2 ../output/alignments/sco_protein_trimmed_alignments/`basename ${fname%%_nt_trimmed.fasta}`_aligned.fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/SCO_backthreaded/`basename ${fname%%_nt_trimmed.fasta}`_backthreaded_aln.fasta
done