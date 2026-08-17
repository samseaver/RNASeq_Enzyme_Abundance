#!/bin/bash
# Rebuild the plastidial reconstructions from the PlantSEED genome-scale models.
#
# The published models were built from ModelSEED/PlantSEED at tag v2.5, whose
# Papers/bioflux-preprint-260807/ directory carries one <Species>_model.json per
# species. Point PLANTSEED at a checkout of that tag:
#
#     git clone https://github.com/ModelSEED/PlantSEED.git
#     git -C PlantSEED checkout v2.5
#     PLANTSEED=$PWD/PlantSEED ./example_plastid_extraction.sh
#
# The extraction needs BiochemPy from ModelSEED/ModelSEEDDatabase:
#
#     export PYTHONPATH=/path/to/ModelSEEDDatabase/Libs/Python
#
# This previously ran against Models/*-reconstruction_cleaned.json, which
# predate v2.5 and do not reproduce the published models.
#
# Output filenames are the ones the rest of the pipeline expects. The Poplar
# model is named inconsistently for historical reasons and is left that way,
# because the companion repository refers to it by that name.
set -euo pipefail

: "${PLANTSEED:?set PLANTSEED to a checkout of ModelSEED/PlantSEED at tag v2.5}"
ART="$PLANTSEED/Papers/bioflux-preprint-260807"

if [ ! -d "$ART" ]; then
    echo "ERROR: $ART not found." >&2
    echo "       Is PLANTSEED a PlantSEED checkout, and is it on tag v2.5?" >&2
    exit 1
fi

extract () {  # <source model> <destination>
    if [ ! -f "$1" ]; then
        echo "ERROR: missing $1" >&2
        exit 1
    fi
    ./extract-plastidial-model.py -m "$1" -o "$2"
}

echo '################################################'
echo '# Running extract-plastidial-model.py for Sorghum'
echo '################################################'
extract "$ART/Sbicolor_v3.1.1_model.json" \
        projects/qpsi-plastidial/inputs/Sbicolor-v3.1.1-plastidial-reconstruction.json

echo
echo '###############################################'
echo '# Running extract-plastidial-model.py for Poplar'
echo '###############################################'
extract "$ART/Ptrichocarpa_v4.1_model.json" \
        projects/qpsi-plastidial/inputs/plastidial-Ptrichocarpa-v4.1-reconstruction_fixed.json

echo
echo '####################################################'
echo '# Running extract-plastidial-model.py for Arabidopsis'
echo '####################################################'
extract "$ART/Athaliana_TAIR10_model.json" \
        data/metabolic_models/plastidial_models/Athaliana-plastidial-reconstruction.json
