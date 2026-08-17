# RNASeq_Enzyme_Abundance

Transcript processing and reaction scoring for the plastidial metabolism of
*Sorghum bicolor* (BTx623) and *Populus trichocarpa* (Nisqually-1) under iron
limitation.

## Preprint

**Simulating Iron Deficiency in Plant Plastidial Metabolism With a Flexible
Neural-Mechanistic Hybrid Approach**

<https://doi.org/10.1101/2025.06.10.658179>

That DOI is versioned and still serves the earlier submission; bioRxiv has not
yet posted the revision this repository corresponds to. Until it does, the
title above is the one to cite, and the code and data here match the revised
manuscript rather than the version currently online.

## Context

This repository turns TMM-normalized transcript abundances into per-reaction
enzyme-capacity estimates for a plastidial metabolic reconstruction, and builds
two of the manuscript figures. It is one of three components:

| | |
|---|---|
| Genome annotation and reconstruction | [ModelSEED/PlantSEED](https://github.com/ModelSEED/PlantSEED) at tag **`v2.5`** |
| **Transcript processing and reaction scores** | **this repository** |
| Mechanistic loss function and gradient descent | [samseaver/Biochem_Informed_Mechanistic_Gradient_Opt](https://github.com/samseaver/Biochem_Informed_Mechanistic_Gradient_Opt) |

The plastidial reconstructions in `projects/qpsi-plastidial/inputs/` were built
from the genome-scale models published under `Papers/bioflux-preprint-260807/`
in PlantSEED `v2.5`. Both study species come out at 375 reactions with
`bio1` optimal at 3.339603.

## Two scores

**Reaction score** `r_s` — the absolute enzyme capacity available to a
reaction, computed from the abundances of the transcripts encoding its
subunits. For a complex the score is limited by its least abundant subunit
(sum-min-sum), and the gene that sets it is reported as the `limiting_subunit`.

**Relative reaction score** `r_s-tilde` — `r_s` divided by the total plastidial
protein pool for that condition, so the two species can be compared as
allocation shares rather than raw abundances.

## Environment

Python 3 with `pandas`, `numpy`, `cobra`, `matplotlib`, `seaborn` and `plotly`
(plus `kaleido==0.2.1` for static PNG export from plotly on a headless
machine). The commands below assume a conda/mamba environment named
`bf-runtime`:

```bash
micromamba run -n bf-runtime python <script>
```

`BiochemPy` from
[ModelSEEDDatabase](https://github.com/ModelSEED/ModelSEEDDatabase) is needed
for the plastidial extraction stage only:

```bash
export PYTHONPATH=/path/to/ModelSEEDDatabase/Libs/Python
```

## Pipeline

Each stage writes the inputs of the next. Stages 1 and 2 only need re-running
when the upstream reconstruction changes; stages 3 and 4 when the transcript
data or the models change.

### 1. Reconcile transcript IDs with model gene IDs

Model gene IDs and RNA-seq gene IDs rarely agree on suffixes. This reports the
overlap under a series of regex strippings and, with `-o`, writes a model whose
IDs match the transcript table.

```bash
./model-transcript-mapper.py \
    -m Models/Sbicolor-v3.1.1-reconstruction.json \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv.xz \
    -o Models/Sbicolor-v3.1.1-reconstruction_cleaned.json
```

`./example_transcript_mapping.sh` runs both species.

### 2. Extract the plastidial model

Takes a genome-scale reconstruction and returns the plastid stroma (`_d0`) and
thylakoid (`_y0`) subnetwork, plus the media exchanges and the plastidial
biomass reaction. Which reactions are injected or excluded is declared in
`src/plastidial_model_extraction/parameters.json` — including the twelve
glucosinolate false positives that PlantSEED places in the plastid of species
that do not make glucosinolates.

```bash
git clone https://github.com/ModelSEED/PlantSEED.git
git -C PlantSEED checkout v2.5

PLANTSEED=$PWD/PlantSEED \
PYTHONPATH=/path/to/ModelSEEDDatabase/Libs/Python \
    ./example_plastid_extraction.sh
```

That rebuilds all three plastidial models — Sorghum, Poplar and Arabidopsis —
in place, and reproduces the tracked files byte-for-byte. For a single species:

```bash
PYTHONPATH=/path/to/ModelSEEDDatabase/Libs/Python \
./extract-plastidial-model.py \
    -m PlantSEED/Papers/bioflux-preprint-260807/Sbicolor_v3.1.1_model.json \
    -o projects/qpsi-plastidial/inputs/Sbicolor-v3.1.1-plastidial-reconstruction.json
```

The genome-scale reconstructions under `Models/` predate `v2.5` and will not
reproduce the published plastidial models; they are kept only as the input to
stage 1's example.

### 3. Reaction scores

Reads the models in `projects/qpsi-plastidial/inputs/` and the TMM tables in
`projects/qpsi-plastidial/rnaseq-data/`, and writes
`{Species}_reaction_scores.tsv` to `projects/qpsi-plastidial/integration_results/`.
Species, project and column names come from `parameters.py`; set
`RESULTS_FOLDER` to write somewhere else and diff before overwriting.

```bash
micromamba run -n bf-runtime python generate_reaction_scores.py
```

Output columns: `condition`, `reaction_id`, `reaction_score`,
`limiting_subunit`. Conditions are `Tissue_Treatment_Timepoint`, e.g.
`Leaf_FeLim_7d`.

### 4. Relative reaction scores

Divides each reaction score by the plastidial protein pool for its condition,
using the Arabidopsis plastid-proteome reference and the ortholog mapping.

```bash
./calculate-plastidial-molar-fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv.xz \
    -r projects/qpsi-plastidial/integration_results/Sorghum_reaction_scores.tsv \
    -o projects/qpsi-plastidial/integration_results/Sorghum_reaction_molar_fractions.tsv \
    --orthologs data/orthologs/Ath-Sbi-Orthologs.tsv
```

`./example_molar_fractions.sh` runs both species.

## Figures

Both take no arguments — every default resolves inside this repository — and
each is named after the file it writes.

```bash
micromamba run -n bf-runtime python figures/fig_scatter_rslt.py
micromamba run -n bf-runtime python figures/fig_proteome.py
```

| script | output | figure |
|---|---|---|
| `figures/fig_scatter_rslt.py` | `fig_scatter_rslt.png`, `.html` | Control vs FeLim reaction scores on log axes, 4 rows x 5 timepoints, coloured by I-dist |
| `figures/fig_proteome.py` | `fig_proteome.png` | plastid vs non-plastid abundance densities, with the method-comparison panel |

**I-dist** is the perpendicular distance of a `(Control, FeLim)` pair to the
identity line, measured in log space:
`(log10(FeLim) - log10(Control)) / sqrt(2)`. Because reaction scores are
log-normal over roughly five decades, taking it in log space makes it a scaled
log fold change, independent of how large the two scores are. The top 5% by
`|I-dist|` are outlined in black.

`plotAbundanceDistributions.py` and `plotCombinedMethodComparison.py` provide
the panels that `fig_proteome.py` assembles and can also be run alone. They
share `LOG_FLOOR_QUANTILE` with `fig_scatter_rslt.py` and must stay in step
with it, or the two figures will disagree about which reactions are in the
top 5%.

## Layout

```
parameters.py                     species, project and column configuration
generate_reaction_scores.py       stage 3 driver
extract-plastidial-model.py       stage 2
model-transcript-mapper.py        stage 1
calculate-plastidial-molar-fractions.py   stage 4
example_*.sh                      both-species wrappers for stages 1, 2 and 4
figures/                          the two manuscript figures and their panels
src/plastidial_model_extraction/  plastid subnetwork extraction; parameters.json
src/reaction_scores/              scoring, model and TMM lookup
src/util/                         KBase model classes, protein-weight helpers
projects/qpsi-plastidial/
    inputs/                       plastidial reconstructions (from PlantSEED v2.5)
    rnaseq-data/                  TMM tables
    integration_results/          reaction scores and molar fractions
data/orthologs/                   Arabidopsis-to-species ortholog tables
data/plastid_proteome/            curated Arabidopsis plastid gene list
data/metabolic_models/            media, biomass, and reference reconstructions
Models/                           pre-v2.5 genome-scale reconstructions
```

## Data availability

The plastidial reconstructions for both species are included, as both derive
from published genomes. For Sorghum the repository also carries the
transcript-derived reaction scores and molar fractions for all 11 leaf
conditions, together with the TMM table. The equivalent Poplar files are
available from the corresponding authors on request.
