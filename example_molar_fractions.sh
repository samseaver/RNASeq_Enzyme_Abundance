#!/bin/bash
# Sorghum
./calculate_plastidial_molar_fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv \
    -r projects/qpsi-plastidial/integration_results/Sorghum_reaction_score.tsv \
    -o projects/qpsi-plastidial/integration_results/Sorghum_reaction_molar_fractions.tsv \
   --orthologs data/orthologs/Ath-Sbi-Orthologs.tsv

# Poplar
./calculate_plastidial_molar_fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Poplar_raw_genes_tmm_mean.tsv \
    -r projects/qpsi-plastidial/integration_results/Poplar_reaction_score.tsv \
    -o projects/qpsi-plastidial/integration_results/Poplar_reaction_molar_fractions.tsv \
   --orthologs data/orthologs/Ath-Ptr-Orthologs.tsv
