import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
import sys
import argparse
import types
import pandas as pa
import plotly.graph_objects as go
import plotly.io as pio

from pathlib import Path
project_root = str(Path(__file__).resolve().parent.parent) + '/'
sys.path.append(project_root)
sys.path.append(os.path.join(project_root, 'src', 'reaction_scores'))

import src.reaction_scores.reactionScoresHelper as rsh
from src.reaction_scores.computeMaxReactionScores import (
    find_model_json, find_tmm_file as _find_tmm_model, load_tmm_data,
    SPECIES_SYNONYMS as MODEL_SYNONYMS,
)
from src.reaction_scores.computeScoresAndPredictions import Species

SPECIES_COLOR = {
    'Poplar':  '#2b6cb0',
    'Sorghum': '#c05621',
}
METHOD_STYLE = {
    'max':      dict(dash='solid', symbol='circle',      label=r"$r_{s} \; max$"),
    'sum':      dict(dash='dash',  symbol='square',      label=r"$r_{s} \; sum$"),
    'relative': dict(dash='dot',   symbol='triangle-up', label=r"$\tilde{r}_{s} \; sum$")# label=r"$r_{s}^{r} \; sum$"),
}
# Candidate suffixes tried in order for each method (handles singular/plural naming)
METHOD_FILES = {
    'sum':      ([('_reaction_scores.tsv', 'reaction_score'),
                  ('_reaction_score.tsv',  'reaction_score')]),
    'max':      ([('_reaction_score_max.tsv', 'reaction_score')]),
    'relative': ([('_reaction_molar_fractions.tsv', 'relative_reaction_score')]),
}
FONT = "Helvetica, Arial, sans-serif"
QUANTILE_THRESHOLD = 0.95
# Numeric day values — used for the x-axis so spacing reflects actual time
DAY_NUMERIC = {'2d': 2, '4d': 4, '7d': 7, '14d': 14, '21d': 21}


def resolve_scores_file(data_dir, species_name, method):
    """Return (path, value_col) for the first candidate that exists, or (None, None)."""
    for suffix, value_col in METHOD_FILES[method]:
        path = os.path.join(data_dir, f"{species_name}{suffix}")
        if os.path.exists(path):
            return path, value_col
    return None, None


def compute_scores_on_the_fly(species_name, method, models_dir, tmm_dir, ignore_organellar_roles):
    """Compute reaction scores via reactionScoresHelper (mirrors computeMaxReactionScores.py)."""
    synonyms = MODEL_SYNONYMS.get(species_name, [species_name])
    model_json = find_model_json(models_dir, species_name)
    tmm_file = _find_tmm_model(tmm_dir, species_name)
    species = Species(species_name, synonyms, model_json)
    tmm_df = load_tmm_data(tmm_file, species)
    params = types.SimpleNamespace(ignore_organellar_roles=ignore_organellar_roles)
    csp = types.SimpleNamespace(rnaSeq_id_col='Gene_ID', value_column='value', group_columns=['condition'])
    scores = rsh.compute_model_score(tmm_df, params, species, csp, method=method, verbose=False)
    scores = scores.rename(columns={csp.value_column: 'reaction_score', csp.rnaSeq_id_col: 'limiting_subunit'})
    return scores[['condition', 'reaction_id', 'reaction_score', 'limiting_subunit']]


def load_and_pivot(source, value_col, tissue, treatment, control_id, exclude_days=None):
    """Pivot reaction scores so control and treatment are separate columns.
    `source` may be a file path (str) or a pre-loaded DataFrame.
    Mirrors load_and_pivot() in generate_reaction_scores_figure.py."""
    if isinstance(source, str):
        df = pa.read_csv(source, sep='\t')
    else:
        df = source.copy()
    parts = df['condition'].str.split('_', n=2, expand=True)
    df['tissue'], df['treatment'], df['day'] = parts[0], parts[1], parts[2]
    # Normalize zero-padded day labels so they match the day_order strings
    df['day'] = df['day'].replace({'02d': '2d', '04d': '4d', '07d': '7d'})
    df = df[df['tissue'] == tissue]
    df = df[df['treatment'].isin([treatment, control_id])]
    if exclude_days:
        df = df[~df['day'].isin(exclude_days)]
    pivoted = df.pivot_table(
        index=['reaction_id', 'day'], columns='treatment',
        values=value_col, aggfunc='first',
    ).reset_index()
    return pivoted.dropna(subset=[control_id, treatment])


def compute_top_reaction_counts(data_dir, species_list, method, tissue, treatment, control_id,
                                 exclude_days, day_order, threshold=QUANTILE_THRESHOLD,
                                 models_dir=None, tmm_dir=None, ignore_organellar_roles=None):
    """Replicate the dark-circle highlight logic from generate_reaction_scores_figure.py:
      - Pivot scores so Control / FeLim are columns
      - Normalize by the single global max across ALL species and ALL days
      - Rank absolute distance from the identity line globally
      - Count reactions at or above `threshold` per species per day

    For 'sum' and 'max' methods: if models_dir + tmm_dir are provided, scores are computed
    on-the-fly via reactionScoresHelper (same approach as computeMaxReactionScores.py).
    For 'relative': reads the precomputed _reaction_molar_fractions.tsv file.

    The ranking scope (all species × all days combined, per method) matches
    prepare_data_group() in generate_reaction_scores_figure.py exactly.
    """
    compute_live = method == 'max' and models_dir and tmm_dir

    all_dfs = []
    for sp in species_list:
        if compute_live:
            try:
                print(f"  [compute] {sp} {method} scores via reactionScoresHelper …")
                raw_df = compute_scores_on_the_fly(
                    sp, method, models_dir, tmm_dir, ignore_organellar_roles)
                value_col = 'reaction_score'
                df = load_and_pivot(raw_df, value_col, tissue, treatment, control_id, exclude_days)
            except Exception as exc:
                print(f"  [warn] on-the-fly compute failed for {sp} {method}: {exc}")
                continue
        else:
            path, value_col = resolve_scores_file(data_dir, sp, method)
            if path is None:
                print(f"  [skip] no {method} file found for {sp} in {data_dir}")
                continue
            df = load_and_pivot(path, value_col, tissue, treatment, control_id, exclude_days)
        df['species'] = sp
        all_dfs.append(df)

    if not all_dfs:
        return {sp: pa.DataFrame(columns=['day', 'counts']) for sp in species_list}
    combined = pa.concat(all_dfs, ignore_index=True)

    # Global normalisation — same as generate_reaction_scores_figure.py
    global_max = max(combined[control_id].max(), combined[treatment].max())
    combined['abs_dist'] = (
        (combined[treatment] - combined[control_id]) / global_max
    ).abs()
    # Global rank across all species + all days (pct=True gives 0–1 fraction)
    combined['rxn_dist_quantile'] = combined['abs_dist'].rank(pct=True)

    all_days_df = pa.DataFrame({'day': day_order})

    result = {}
    for sp in species_list:
        sp_df = combined[combined['species'] == sp]
        top_df = sp_df[sp_df['rxn_dist_quantile'] >= threshold]
        counts = top_df.groupby('day').size().reset_index(name='counts')
        # Left-join so days with zero reactions appear explicitly as 0
        counts = all_days_df.merge(counts, on='day', how='left').fillna(0)
        counts['counts'] = counts['counts'].astype(int)
        counts['day'] = pa.Categorical(counts['day'], categories=day_order, ordered=True)
        result[sp] = counts.sort_values('day')

    return result


def add_species_traces(fig, species_name, method_counts, day_order):
    """Add one line+marker trace per method for the given species."""
    color = SPECIES_COLOR[species_name]
    max_count = 0

    for method in ['max', 'sum', 'relative']:
        counts_df = method_counts[method].get(species_name, pa.DataFrame())
        if not counts_df.empty:
            max_count = max(max_count, counts_df['counts'].max())

        m = METHOD_STYLE[method]
        x_vals = counts_df['day'].map(DAY_NUMERIC) if not counts_df.empty else []
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=counts_df['counts'] if not counts_df.empty else [],
            mode='lines+markers',
            line=dict(color=color, dash=m['dash'], width=2.2),
            marker=dict(symbol=m['symbol'], size=9, color='white',
                        line=dict(color=color, width=2.2)),
            showlegend=False,
        ))

    return max_count


def add_legend_traces(fig, species_list):
    """Split legend: colored solid lines for species identity,
    then black lines+markers for scoring method."""
    for species_name in species_list:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            name=species_name,
            mode='lines',
            line=dict(color=SPECIES_COLOR[species_name], width=2.5, dash='solid'),
            showlegend=True,
        ))
    for method in ['max', 'sum', 'relative']:
        m = METHOD_STYLE[method]
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            name=m['label'],
            mode='lines+markers',
            line=dict(color='black', dash=m['dash'], width=2),
            marker=dict(symbol=m['symbol'], size=8, color='white',
                        line=dict(color='black', width=2)),
            showlegend=True,
        ))


def plot_combined_figure(data_dir, species_list, treatment='FeLim', control_id='Control', tissue='Leaf',
                         day_order=('2d', '4d', '7d', '14d', '21d'), exclude_days=None,
                         output_path='combined_method_comparison.png', scale=4,
                         models_dir=None, tmm_dir=None, ignore_organellar_roles=None):
    # Pre-compute counts for every method using the global ranking approach
    method_counts = {}
    for method in ['max', 'sum', 'relative']:
        method_counts[method] = compute_top_reaction_counts(
            data_dir, species_list, method, tissue, treatment, control_id,
            exclude_days, list(day_order),
            models_dir=models_dir, tmm_dir=tmm_dir,
            ignore_organellar_roles=ignore_organellar_roles,
        )

    fig = go.Figure()

    max_count = 0
    for species_name in species_list:
        panel_max = add_species_traces(fig, species_name, method_counts, list(day_order))
        max_count = max(max_count, panel_max)

    add_legend_traces(fig, species_list)

    fig.update_yaxes(range=[-2, max(max_count * 1.12, 1)])

    fig.update_layout(
        template='plotly_white',
        paper_bgcolor='white',
        plot_bgcolor='white',
        height=300,
        width=520,
        font=dict(family=FONT, size=13),
        legend=dict(
            orientation="v",
            yanchor="top", y=0.69,
            xanchor="right", x=0.47,
            bgcolor='rgba(255,255,255,0.85)',
            bordercolor='rgba(0,0,0,0.15)',
            borderwidth=1,
            font=dict(size=9),
            tracegroupgap=0,
        ),
        margin=dict(l=60, r=15, t=15, b=45),
    )
    day_vals = [DAY_NUMERIC[d] for d in day_order if d in DAY_NUMERIC]
    fig.update_xaxes(
        tickmode='array',
        tickvals=day_vals,
        ticktext=[d for d in day_order if d in DAY_NUMERIC],
        tickangle=0,
        tickfont=dict(size=12),
        showgrid=False,
        linecolor='rgba(0,0,0,0.3)',
        linewidth=1,
    )
    fig.update_yaxes(
        tickfont=dict(size=12),
        gridcolor='rgba(0,0,0,0.1)',
        griddash='dot',
        title_text="# reactions at 95th percentile",
        title_font=dict(size=11),
        title_standoff=6,
        zeroline=False,
    )

    output_format = os.path.splitext(output_path)[1].lstrip('.').lower() or 'png'
    img_bytes = pio.to_image(fig, format=output_format, scale=scale)
    with open(output_path, "wb") as f:
        f.write(img_bytes)

    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a single Poplar+Sorghum reaction-score method comparison figure from precomputed TSVs."
    )
    parser.add_argument("--data-dir",
                         default=os.path.join(project_root, 'projects',
                                               'qpsi-plastidial', 'integration_results'))
    parser.add_argument("--species", nargs='+', default=["Poplar", "Sorghum"])
    parser.add_argument("--treatment", default="FeLim")
    parser.add_argument("--control", default="Control")
    parser.add_argument("--tissue", default="Leaf")
    parser.add_argument("--exclude-days", nargs='*', default=["0h", "1h"])
    parser.add_argument("--output", default="combined_method_comparison.png")
    parser.add_argument("--scale", type=float, default=4)
    parser.add_argument("--models-dir",
                         default=os.path.join(project_root, "Models"),
                         help="Folder with plastidial-reconstruction model JSON files")
    parser.add_argument("--tmm-dir",
                         default=os.path.join(project_root, "data", "RNAseq", "tmm"),
                         help="Folder with TMM CSV files (used to compute sum/max scores on-the-fly)")
    parser.add_argument("--ignore-organellar-roles",
                         default=os.path.join(project_root, "data",
                                               "organellar-encoded_subunits_to_ignore.txt"))
    args = parser.parse_args()

    plot_combined_figure(args.data_dir, args.species, args.treatment, args.control, args.tissue,
                         exclude_days=args.exclude_days, output_path=args.output, scale=args.scale,
                         models_dir=args.models_dir, tmm_dir=args.tmm_dir,
                         ignore_organellar_roles=args.ignore_organellar_roles)
