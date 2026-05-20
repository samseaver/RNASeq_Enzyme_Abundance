#!/bin/bash
echo '#############################################################'
echo '# Running calculate-plastidial-molar-fractions.py for Sorghum'
echo '#############################################################'
./calculate-plastidial-molar-fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv \
    -r projects/qpsi-plastidial/integration_results/Sorghum_reaction_score.tsv \
    -o projects/qpsi-plastidial/integration_results/Sorghum_reaction_molar_fractions.tsv \
   --orthologs data/orthologs/Ath-Sbi-Orthologs.tsv
echo
echo
echo '############################################################'
echo '# Running calculate-plastidial-molar-fractions.py for Poplar'
echo '############################################################'
./calculate-plastidial-molar-fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Poplar_raw_genes_tmm_mean.tsv \
    -r projects/qpsi-plastidial/integration_results/Poplar_reaction_score.tsv \
    -o projects/qpsi-plastidial/integration_results/Poplar_reaction_molar_fractions.tsv \
   --orthologs data/orthologs/Ath-Ptr-Orthologs.tsv
