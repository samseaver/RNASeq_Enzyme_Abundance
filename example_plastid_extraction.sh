#!/bin/bash
echo '################################################'
echo '# Running extract-plastidial-model.py for Sorghum'
echo '################################################'
./extract-plastidial-model.py \
    -m Models/Sbicolor-v3.1.1-reconstruction_cleaned.json \
    -o Models/Sbicolor-v3.1.1-plastidial-reconstruction.json

echo
echo
echo '###############################################'
echo '# Running extract-plastidial-model.py for Poplar'
echo '###############################################'
./extract-plastidial-model.py \
    -m Models/Ptrichocarpa-v4.1-reconstruction_cleaned.json \
    -o Models/Ptrichocarpa-v4.1-plastidial-reconstruction.json
