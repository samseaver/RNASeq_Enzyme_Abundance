#!/bin/bash
echo '################################################'
echo '# Running model-transcript-mapper.py for Sorghum'
echo '################################################'
./model-transcript-mapper.py \
    -m Models/Sbicolor-v3.1.1-reconstruction.json \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv.xz \
    -o Models/Sbicolor-v3.1.1-reconstruction_cleaned.json \
    -r '\.\d+\.p$'
echo
echo
echo '###############################################'
echo '# Running model-transcript-mapper.py for Poplar'
echo '###############################################'
./model-transcript-mapper.py \
    -m Models/Ptrichocarpa-v4.1-reconstruction.json \
    -t projects/qpsi-plastidial/rnaseq-data/Poplar_raw_genes_tmm_mean.tsv \
    -o Models/Ptrichocarpa-v4.1-reconstruction_cleaned.json \
    -r '\.\d+\.p$'
