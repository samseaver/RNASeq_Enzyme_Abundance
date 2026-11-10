import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
import sys
import argparse
import types
import pandas as pa

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)

import src.reaction_scores.reactionScoresHelper as rsh
from src.reaction_scores.computeScoresAndPredictions import Species

# Same species -> filename-synonym mapping as parameters.py, used to locate
# each species' plastidial-reconstruction model JSON under --models-dir.
SPECIES_SYNONYMS = {
    "Poplar": ["Poplar", "Pptrich", "Ptrichocarpa", "ptrich", "ptr"],
    "Sorghum": ["Sorghum", "Sbicolor", "Sbi", "sbi"],
}



def find_model_json(models_dir, species_name):
    # Match "plastidial" and "reconstruction" independently rather than as the
    # single substring "plastidial-reconstruction". The two models are named
    # inconsistently -- Sbicolor-v3.1.1-plastidial-reconstruction.json but
    # plastidial-Ptrichocarpa-v4.1-reconstruction_fixed.json -- and requiring
    # the joined form silently dropped Poplar.
    synonyms = SPECIES_SYNONYMS.get(species_name, [species_name])
    candidates = []
    for file_name in os.listdir(models_dir):
        if not file_name.endswith('.json'):
            continue
        lower = file_name.lower()
        if "plastidial" not in lower or "reconstruction" not in lower:
            continue
        if "cleaned" in lower or "media" in lower:
            continue
        if any(syn.lower() in lower for syn in synonyms):
            candidates.append(file_name)
    if not candidates:
        raise FileNotFoundError(
            f"No plastidial reconstruction JSON found for {species_name} in {models_dir}")
    if len(candidates) > 1:
        # Silently taking sorted()[0] is how a stale model gets picked up
        # without anyone noticing. Make the caller choose.
        raise ValueError(
            f"Ambiguous plastidial reconstruction for {species_name} in {models_dir}: "
            f"{sorted(candidates)}. Pass --models-dir to a directory with one per species.")
    return os.path.join(models_dir, candidates[0])


def find_tmm_file(tmm_dir, species_name):
    synonyms = SPECIES_SYNONYMS.get(species_name, [species_name])
    for file_name in os.listdir(tmm_dir):
        if any(syn.lower() in file_name.lower() for syn in synonyms):
            return os.path.join(tmm_dir, file_name)
    raise FileNotFoundError(f"No TMM file found for {species_name} in {tmm_dir}")


# Mirrors computeScoresAndPredictions.readTMMdata: restrict to genes present in
# the model. TMM outlier capping removed --- values are used uncapped.
def load_tmm_data(tmm_file, species, id_col='Gene_ID', value_col='value'):
    # CLAUDE 2026-08-12: the TMM tables under projects/*/rnaseq-data are
    # tab-separated and already carry a `condition` column of the form
    # Tissue_Treatment_Timepoint. Sniff the delimiter, and only rebuild
    # `condition` when the three source columns are actually present.
    if tmm_file.endswith(('.xz', '.gz', '.bz2')):
        sep = '\t'
    else:
        with open(tmm_file, 'rb') as fh:
            sep = '\t' if fh.readline().count(b'\t') else ','
    tmm_df = pa.read_csv(tmm_file, sep=sep)

    if tmm_df.columns[0] != id_col:
        tmm_df.columns.values[0] = id_col

    if 'condition' not in tmm_df.columns:
        tmm_df['condition'] = (tmm_df['tissue'].astype(str) + '_'
                                + tmm_df['treatment'].astype(str) + '_'
                                + tmm_df['time_stamp'].astype(str))

    tmm_df = tmm_df[tmm_df[id_col].isin(species.metModel.modelfeatures_dict)]

    # TMM outlier capping removed --- values are used uncapped.

    return tmm_df


def compute_max_reaction_scores(species_name, models_dir, tmm_dir, ignore_organellar_roles, verbose=False):
    model_json = find_model_json(models_dir, species_name)
    tmm_file = find_tmm_file(tmm_dir, species_name)

    if verbose:
        print(f"{species_name}: model={model_json}")
        print(f"{species_name}: tmm={tmm_file}")

    species = Species(species_name, SPECIES_SYNONYMS.get(species_name, [species_name]), model_json)

    tmm_df = load_tmm_data(tmm_file, species)

    params = types.SimpleNamespace(ignore_organellar_roles=ignore_organellar_roles)
    csp = types.SimpleNamespace(rnaSeq_id_col='Gene_ID', value_column='value', group_columns=['condition'])

    reaction_scores = rsh.compute_model_score(tmm_df, params, species, csp, method='max', verbose=verbose)

    reaction_scores = reaction_scores.rename(columns={
        csp.value_column: 'reaction_score',
        csp.rnaSeq_id_col: 'limiting_subunit',
    })
    reaction_scores = reaction_scores[['condition', 'reaction_id', 'reaction_score', 'limiting_subunit']]

    return reaction_scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute r_s max reaction scores (method='max') via reactionScoresHelper.compute_model_score, "
                     "matching the approach used in computeScoresAndPredictions.generate_reactionScores."
    )
    parser.add_argument("--species", nargs='+', default=["Poplar", "Sorghum"])
    parser.add_argument("--models-dir",
                         default=os.path.join(project_root, "projects", "qpsi-plastidial", "inputs"))
    parser.add_argument("--tmm-dir",
                         default=os.path.join(project_root, "projects", "qpsi-plastidial", "rnaseq-data"))
    parser.add_argument("--ignore-organellar-roles",
                         default=os.path.join(project_root, "data", "organellar-encoded_subunits_to_ignore.txt"))
    parser.add_argument("--output-dir", default=os.path.join(project_root, "data", "other_input_files"))
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for species_name in args.species:
        reaction_scores = compute_max_reaction_scores(
            species_name, args.models_dir, args.tmm_dir, args.ignore_organellar_roles, verbose=args.verbose)

        out_path = os.path.join(args.output_dir, f"{species_name}_reaction_score_max.tsv")
        reaction_scores.to_csv(out_path, sep='\t', index=False)
        print(f"Saved {species_name} r_s max scores to {out_path}")
