# RNASeq_Enzyme_Abundance

Transcript processing and reaction scoring for the plastidial metabolism of
*Sorghum bicolor* (BTx623) and *Populus trichocarpa* (Nisqually-1) under iron
limitation.

## Preprint

**Simulating Iron Deficiency in Plant Plastidial Metabolism With a Flexible
Neural-Mechanistic Hybrid Approach**

El Alaoui, S., Henry, C. S., Blaby-Haas, C., Paape, T., Xie, M., and
Seaver, S. M. bioRxiv, version 2, posted 2026-08-17.
<https://doi.org/10.1101/2025.06.10.658179>

## Context

This repository turns TMM-normalized transcript abundances into per-reaction
enzyme-capacity estimates for a plastidial metabolic reconstruction, and builds
two figures for the preprint. It is one of three components:

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

**Reaction score** $r_s$ — the absolute enzyme capacity available to a
reaction, computed from the abundances of the transcripts encoding its
subunits. For a complex the score is limited by its least abundant subunit
(sum-min-sum), and the gene that sets it is reported as the `limiting_subunit`.

**Relative reaction score** $\tilde{r}_s$ — $r_s$ divided by the total
plastidial protein pool for that condition, so the two species can be compared
as allocation shares rather than raw abundances.

*Why the second score exists.* The two species do not put comparable amounts of
transcript into the plastid to begin with, so their absolute scores are not on a
common footing. Summed over plastid-localized genes in leaf
(`integration_results/<Species>_plastid_transcript_totals.tsv`):

| | mean | range across leaf conditions |
|---|---|---|
| Sorghum | 836,000 | 327,000 – 1,259,000 |
| Poplar | 495,000 | 357,000 – 611,000 |

Sorghum's plastid pool is not only larger on average but swings roughly
four-fold across the time course, against under two-fold for Poplar. Comparing
raw $r_s$ between them, or across timepoints within Sorghum, largely measures
that pool rather than the enzyme. The normalization makes the comparison
meaningful, and it is the score behind the cross-species allocation results.

*What it must not be used for.* $\tilde{r}_s$ is a **share of a moving total**,
not a capacity. A reaction's share can rise while its absolute enzyme abundance
falls, simply because the pool around it contracted faster — which is exactly
what happens in Sorghum after day 7. It therefore cannot be used to bound or
interpret flux. The flux work uses the un-normalized $r_s$, and the companion
repository reads `<Species>_reaction_scores.tsv` and not the molar-fractions
file for precisely this reason. Use $\tilde{r}_s$ for allocation, $r_s$ for
capacity.

## Environment

Python 3 with `pandas`, `numpy`, `cobra`, `matplotlib`, `seaborn` and `plotly`
The commands below assume a micromamba environment named `bf-runtime`:

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

Three stages. Each writes the inputs of the next.

| stage | reads | writes |
|---|---|---|
| 1. Extract the plastidial model | a PlantSEED `v2.5` genome-scale model | `projects/qpsi-plastidial/inputs/<Species>-…-reconstruction.json` |
| 2. Reaction scores | that model + `rnaseq-data/<Species>_raw_genes_tmm_mean.tsv[.xz]` | `integration_results/<Species>_reaction_scores.tsv` |
| 3. Relative reaction scores | those scores + the TMM table + the plastid proteome reference | `integration_results/<Species>_reaction_molar_fractions.tsv` |

Stage 1 only needs re-running when the upstream reconstruction changes; stages
2 and 3 when the transcript data or the models change. Nothing in the pipeline
depends on the preliminary checks described further down.

### 1. Extract the plastidial model

Takes a genome-scale reconstruction and returns the plastid stroma (`_d0`) and
thylakoid (`_y0`) subnetwork, plus the media exchanges and the plastidial
biomass reaction. Which reactions are injected or excluded is declared in
`src/plastidial_model_extraction/parameters.json`.

**This stage rests on Arabidopsis orthology that is already in place.** Both the
plastid and thylakoid localizations it selects on, and the gene-protein-reaction
rules it carries through, come from PlantSEED's OrthoFinder-based annotation,
which propagates curated Arabidopsis roles and compartments onto the target
genome. All 375 reactions arrive with their GPRs already attached. The
extraction itself only works on the individual models that already have GPRs —
that work is finished in the `v2.5` artifacts, against the reference set of
species those artifacts were built for.

Applying this work to a new species requires PlantSEED annotation of that
genome first, which means OrthoFinder orthologs against the curated Arabidopsis
reference.

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

### 2. Reaction scores

Reads the models in `projects/qpsi-plastidial/inputs/` and the TMM tables in
`projects/qpsi-plastidial/rnaseq-data/`, and writes
`{Species}_reaction_scores.tsv` to `projects/qpsi-plastidial/integration_results/`.
Species, project and column names come from `parameters.py`.

```bash
micromamba run -n bf-runtime python generate_reaction_scores.py
```

Output columns: `condition`, `reaction_id`, `reaction_score`,
`limiting_subunit`. Conditions are `Tissue_Treatment_Timepoint`, e.g.
`Leaf_FeLim_7d`.

### 3. Relative reaction scores

Divides each reaction score by the plastidial protein pool for its condition. A
species gene counts towards that pool if any of its Arabidopsis orthologs is in
the curated plastid-proteome list, so this stage needs an Arabidopsis ortholog
table for the species.

Those tables ship with the repository, one per study species, and come from the
same OrthoFinder run as the annotation:

    data/orthologs/Ath-Sbi-Orthologs.tsv
    data/orthologs/Ath-Ptr-Orthologs.tsv

They are also what `figures/fig_proteome.py` reads, through `--orthologs-dir`. A
new species needs its own equivalent table against Arabidopsis.

The output of this stage is for comparing allocation between species and across
time, not for constraining flux — see "Two scores" above.

```bash
./calculate-plastidial-molar-fractions.py \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv.xz \
    -r projects/qpsi-plastidial/integration_results/Sorghum_reaction_scores.tsv \
    -o projects/qpsi-plastidial/integration_results/Sorghum_reaction_molar_fractions.tsv \
    --orthologs data/orthologs/Ath-Sbi-Orthologs.tsv
```

`./example_molar_fractions.sh` runs both species.

## Preliminary check: do the transcript IDs match the model?

Not a pipeline stage. Reaction scoring silently produces nothing for a gene
whose ID in the model does not match its ID in the transcript table, and gene
IDs routinely differ by a transcript suffix or an assembly-version prefix. This
reports the overlap before you spend time on a run that would come back empty.

```bash
micromamba run -n bf-runtime python model-transcript-mapper.py \
    -m projects/qpsi-plastidial/inputs/Sbicolor-v3.1.1-plastidial-reconstruction.json \
    -t projects/qpsi-plastidial/rnaseq-data/Sorghum_raw_genes_tmm_mean.tsv.xz
```

On the tracked Sorghum model that reports 398 of 416 model genes matched. The
18 that do not match are genes with no transcript in the dataset, not an ID
mismatch — which is the distinction the report exists to make. It tries a series
of regex strippings to find a systematic suffix difference, and says so when
there isn't one. Passing `-o` writes a model with the ID fixes applied.

`test_model_with_fba.py` is the other check in the same spirit: it loads a model
in KBase or COBRA format and runs FBA/FVA, to confirm the network still solves.
It needs `cobrakbase` for the KBase format.

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

**I-dist** is the perpendicular distance of a (Control, FeLim) pair to the
identity line, measured in log space:

$$\text{I-dist} = \frac{\log_{10}(\text{FeLim}) - \log_{10}(\text{Control})}{\sqrt{2}}$$

Because reaction scores are log-normal over roughly five decades, taking the
distance in log space makes I-dist a scaled log fold change, independent of how
large the two scores are. The top 5% by $|\text{I-dist}|$ are outlined in
black.

`plotAbundanceDistributions.py` and `plotCombinedMethodComparison.py` provide
the panels that `fig_proteome.py` assembles and can also be run alone. They
share `LOG_FLOOR_QUANTILE` with `fig_scatter_rslt.py` and must stay in step
with it, or the two figures will disagree about which reactions are in the
top 5%.

## Layout

```
parameters.py                     species, project and column configuration
extract-plastidial-model.py       stage 1
generate_reaction_scores.py       stage 2 driver
calculate-plastidial-molar-fractions.py   stage 3
example_*.sh                      both-species wrappers for stages 1 and 3
model-transcript-mapper.py        preliminary ID check (see above)
test_model_with_fba.py            loads a model and runs FBA/FVA on it
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
```

## Data availability

The plastidial reconstructions for both species are included, as both derive
from published genomes. For Sorghum the repository also carries the
transcript-derived reaction scores and molar fractions for all 11 leaf
conditions, together with the TMM table. The equivalent Poplar files are
available from the corresponding authors on request.
