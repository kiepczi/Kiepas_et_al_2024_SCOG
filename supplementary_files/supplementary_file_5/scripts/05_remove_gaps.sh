#!/bin/bash
# trim_gaps.sh
# This scrip was used to trim the alignmets gaps with trimal with -colnumbering parameter and to generate file with the matchin columns in the new and old alignment


trimal -in ../output/alignments/concatenated/initial_concatenated_sco.fasta -out ../output/alignments/concatenated/no_gaps_concatenated_sco.fasta -automated1 -colnumbering > ../output/alignments/concatenated/no_gaps_concatenated_columns.txt