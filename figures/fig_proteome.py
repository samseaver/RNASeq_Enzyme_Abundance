"""
Figure 3: plastid proteome distributions paired with the method comparison.

  Panel A  –  per-day KDE density grid (species x days), plastid vs
              non-plastid genes, Control solid / FeLim dashed
  Panel B  –  reaction-score method comparison line plot: # scored reactions
              in the 95th percentile of |I-dist| over time,
              blue = Poplar / orange = Sorghum

Imports data-loading functions directly from the two source scripts so all
logic stays in one place and this script stays thin.

Data sources
------------
projects/qpsi-plastidial/integration_results/{Species}_reaction_score*.tsv
projects/qpsi-plastidial/rnaseq-data/{Species}_raw_genes_tmm_mean.tsv[.xz]
data/plastid_proteome/compiled_arabidopsis_plastid_genes.txt
a directory of plastidial reconstruction JSONs (--models-dir)

Outputs
-------
fig_proteome.png (override with --output).

Two defaults do not match this repository's layout and must be passed:

    cd figures
    micromamba run -n bf-runtime python fig_proteome.py \
      --tmm-dir ../projects/qpsi-plastidial/rnaseq-data \
      --models-dir <dir of *plastidial-reconstruction*.json>

--tmm-dir defaults to data/RNAseq/tmm, which does not exist here. --models-dir
cannot yet point at projects/qpsi-plastidial/inputs: find_model_json matches the
literal substring "plastidial-reconstruction", and the Poplar model is named
plastidial-Ptrichocarpa-v4.1-reconstruction_fixed.json, so Poplar is skipped.
"""
import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
import sys
import argparse
import pandas as pa
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import seaborn as sns

from pathlib import Path
project_root = str(Path(__file__).resolve().parent.parent) + '/'
sys.path.append(project_root)
sys.path.append(str(Path(__file__).resolve().parent))

from plotAbundanceDistributions import (
    load_classified_tmm, PLASTID_COLOR, NON_PLASTID_COLOR,
)
from plotCombinedMethodComparison import (
    compute_top_reaction_counts, SPECIES_COLOR, METHOD_STYLE, DAY_NUMERIC,
)

FONT       = "Helvetica Neue, Helvetica, Arial, sans-serif"
DAY_ORDER  = ('2d', '4d', '7d', '14d', '21d')
LS_MAP     = {'solid': '-', 'dash': '--', 'dot': ':'}
MARKER_MAP = {'circle': 'o', 'square': 's', 'triangle-up': '^'}


# ── Panel A helpers ──────────────────────────────────────────────────────────

def _draw_density_grid(axes_grid, species_data, species_list, day_order):
    days = list(day_order)
    n_rows = len(species_list)

    for row_idx, species_name in enumerate(species_list):
        df = species_data[species_name]
        for col_idx, day in enumerate(days):
            ax = axes_grid[row_idx][col_idx]
            day_df = df[df['time_stamp'] == day]

            for treatment, ls in [('Control', '-'), ('FeLim', '--')]:
                trmt_df = day_df[day_df['treatment'] == treatment]
                if trmt_df.empty:
                    continue
                plastid_df    = trmt_df[trmt_df['is_plastid']]
                nonplastid_df = trmt_df[~trmt_df['is_plastid']]
                if len(plastid_df) > 5:
                    sns.kdeplot(data=plastid_df, x='value', log_scale=True,
                                color=PLASTID_COLOR, linestyle=ls,
                                linewidth=1.4, ax=ax)
                if len(nonplastid_df) > 5:
                    sns.kdeplot(data=nonplastid_df, x='value', log_scale=True,
                                color=NON_PLASTID_COLOR, linestyle=ls,
                                linewidth=1.4, ax=ax)

            ax.set_facecolor('white')
            ax.tick_params(labelsize=7)

            # Column titles (days) — top row only
            if row_idx == 0:
                ax.set_title(day, fontsize=9, fontweight='bold', pad=3)

            # Row labels (species) — leftmost column only
            if col_idx == 0:
                ax.set_ylabel(species_name, fontsize=9, fontweight='bold', labelpad=4)
            else:
                ax.set_ylabel('')

            # x-axis label — bottom row only
            if row_idx == n_rows - 1:
                ax.set_xlabel('TMM abundance', fontsize=8)
            else:
                ax.set_xlabel('')

            # Suppress inner tick labels (shared axes handle range; hide clutter)
            if col_idx != 0:
                ax.tick_params(labelleft=False)
            if row_idx != n_rows - 1:
                ax.tick_params(labelbottom=False)

            ax.set_ylim(0, 0.7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)


def _density_legend(fig, bbox):
    """Return a figure-level legend for the density panel."""
    handles = [
        Line2D([0], [0], color=PLASTID_COLOR,    lw=1.5, ls='-',  label='Plastid (Control)'),
        Line2D([0], [0], color=PLASTID_COLOR,    lw=1.5, ls='--', label='Plastid (FeLim)'),
        Line2D([0], [0], color=NON_PLASTID_COLOR, lw=1.5, ls='-',  label='Non-plastid (Control)'),
        Line2D([0], [0], color=NON_PLASTID_COLOR, lw=1.5, ls='--', label='Non-plastid (FeLim)'),
    ]
    return fig.legend(handles=handles, loc='upper right', bbox_to_anchor=bbox,
                      fontsize=8, frameon=True, framealpha=0.9,
                      labelspacing=0.2, handletextpad=0.4,
                      handlelength=1.6, borderpad=0.35, ncol=2)


# ── Panel B helpers ──────────────────────────────────────────────────────────

def _draw_line_panel(ax, method_counts, species_list, day_order):
    days     = list(day_order)
    x_ticks  = [DAY_NUMERIC[d] for d in days if d in DAY_NUMERIC]

    for species_name in species_list:
        color = SPECIES_COLOR[species_name]
        for method in ['max', 'sum', 'relative']:
            m          = METHOD_STYLE[method]
            counts_df  = method_counts[method].get(species_name, pa.DataFrame())
            if counts_df.empty:
                continue
            x = counts_df['day'].map(DAY_NUMERIC).values
            y = counts_df['counts'].values
            ax.plot(x, y,
                    color=color,
                    linestyle=LS_MAP[m['dash']],
                    linewidth=1.8,
                    marker=MARKER_MAP[m['symbol']],
                    markersize=7,
                    markerfacecolor='white',
                    markeredgecolor=color,
                    markeredgewidth=1.8)

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(days, fontsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xlim(x_ticks[0] - 1, x_ticks[-1] + 1)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('# reactions at 95th percentile', fontsize=9)
    ax.set_xlabel('Time point', fontsize=9)
    ax.set_axisbelow(True)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _line_legend(ax, species_list):
    handles = []
    # Species — colored solid lines
    for sp in species_list:
        handles.append(Line2D([0], [0], color=SPECIES_COLOR[sp],
                               lw=2, label=sp))
    # Methods — black lines with marker + style
    for method in ['max', 'sum', 'relative']:
        m = METHOD_STYLE[method]
        handles.append(Line2D(
            [0], [0], color='black',
            lw=1.5, linestyle=LS_MAP[m['dash']],
            marker=MARKER_MAP[m['symbol']],
            markerfacecolor='white', markeredgecolor='black',
            markersize=7, label=m['label'],
        ))
    # Tucked against the y-axis so the legend fills the gap between it and the
    # 7d peak. The x-limits run day 1 to day 22, so 7d sits at (7 - 1) / 21 =
    # 0.286 of the axis; the legend is ~0.19 wide, so anchoring its left edge
    # at 0.01 puts its centre near 0.105 and its right edge clear of the peak.
    ax.legend(handles=handles, fontsize=8, loc='upper left',
              frameon=True, framealpha=0.9,
              labelspacing=0.2, handletextpad=0.4,
              handlelength=1.8, borderpad=0.35,
              bbox_to_anchor=(0.01, 0.99))


# ── Main assembly ─────────────────────────────────────────────────────────────

def build_combined_figure(
    data_dir, tmm_dir, orthologs_dir, plastid_gene_file,
    species_list=('Poplar', 'Sorghum'),
    treatment='FeLim', control_id='Control', tissue='Leaf',
    day_order=DAY_ORDER, exclude_days=('0h', '1h'),
    output_path='fig_proteome.png',
    models_dir=None, ignore_organellar_roles=None,
):
    # ── Load data ─────────────────────────────────────────────────────────
    print("Loading abundance data for panel A …")
    species_data = {}
    for sp in species_list:
        species_data[sp] = load_classified_tmm(
            sp, tmm_dir, orthologs_dir, plastid_gene_file,
            tissue=tissue, treatments=(treatment, control_id),
            exclude_days=exclude_days,
        )

    print("Computing reaction-score counts for panel B …")
    method_counts = {}
    for method in ['max', 'sum', 'relative']:
        method_counts[method] = compute_top_reaction_counts(
            data_dir, list(species_list), method, tissue, treatment, control_id,
            list(exclude_days), list(day_order),
            models_dir=models_dir, tmm_dir=tmm_dir,
            ignore_organellar_roles=ignore_organellar_roles,
        )

    # ── Build figure ──────────────────────────────────────────────────────
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial']

    n_days    = len(day_order)
    n_species = len(species_list)

    fig = plt.figure(figsize=(13, 8), facecolor='white')

    # Outer grid: Panel A (top, taller) and Panel B (bottom)
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[n_species * 1.6, 1],
        hspace=0.17,
    )

    # Panel A — density subgrid
    inner_top = outer[0].subgridspec(n_species, n_days, hspace=0.15, wspace=0.08)
    axes_top  = [[fig.add_subplot(inner_top[r, c]) for c in range(n_days)]
                  for r in range(n_species)]

    # Share x and y across ALL density panels
    ax_ref = axes_top[0][0]
    for r in range(n_species):
        for c in range(n_days):
            if not (r == 0 and c == 0):
                axes_top[r][c].sharex(ax_ref)
                axes_top[r][c].sharey(ax_ref)

    _draw_density_grid(axes_top, species_data, list(species_list), day_order)
    _density_legend(fig, bbox=(0.4, 0.6))

    # Panel B — line plot
    ax_line = fig.add_subplot(outer[1])
    _draw_line_panel(ax_line, method_counts, list(species_list), day_order)
    _line_legend(ax_line, list(species_list))

    # Panel labels
    axes_top[0][0].text(-0.18, 1.18, 'A', transform=axes_top[0][0].transAxes,
                         fontsize=14, fontweight='bold', va='top')
    ax_line.text(-0.07, 1.12, 'B', transform=ax_line.transAxes,
                  fontsize=14, fontweight='bold', va='top')

    # ── Save ──────────────────────────────────────────────────────────────
    fmt = os.path.splitext(output_path)[1].lstrip('.').lower() or 'pdf'
    fig.savefig(output_path, dpi=300, facecolor='white',
                bbox_inches='tight', format=fmt)
    plt.close(fig)
    print(f"Saved {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Build combined panel figure (density grid + method line plot).'
    )
    parser.add_argument('--data-dir',
                         default=os.path.join(project_root, 'projects',
                                               'qpsi-plastidial', 'integration_results'),
                         help='Folder with {Species}_reaction_score*.tsv files')
    parser.add_argument('--tmm-dir',
                         default=os.path.join(project_root, 'data', 'RNAseq', 'tmm'))
    parser.add_argument('--orthologs-dir',
                         default=os.path.join(project_root, 'data', 'orthologs'))
    parser.add_argument('--plastid-gene-file',
                         default=os.path.join(project_root, 'data', 'plastid_proteome',
                                               'compiled_arabidopsis_plastid_genes.txt'))
    parser.add_argument('--species', nargs='+', default=['Poplar', 'Sorghum'])
    parser.add_argument('--treatment', default='FeLim')
    parser.add_argument('--control',   default='Control')
    parser.add_argument('--tissue',    default='Leaf')
    parser.add_argument('--exclude-days', nargs='*', default=['0h', '1h'])
    parser.add_argument('--output', default='fig_proteome.png')
    parser.add_argument('--models-dir',
                         default=os.path.join(project_root, 'Models'),
                         help='Folder with plastidial-reconstruction model JSON files')
    parser.add_argument('--ignore-organellar-roles',
                         default=os.path.join(project_root, 'data',
                                               'organellar-encoded_subunits_to_ignore.txt'))
    args = parser.parse_args()

    build_combined_figure(
        data_dir        = args.data_dir,
        tmm_dir         = args.tmm_dir,
        orthologs_dir   = args.orthologs_dir,
        plastid_gene_file = args.plastid_gene_file,
        species_list    = args.species,
        treatment       = args.treatment,
        control_id      = args.control,
        tissue          = args.tissue,
        exclude_days    = tuple(args.exclude_days),
        output_path     = args.output,
        models_dir      = args.models_dir,
        ignore_organellar_roles = args.ignore_organellar_roles,
    )
