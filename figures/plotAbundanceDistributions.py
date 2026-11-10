import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
import sys
import argparse
import pandas as pa
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
project_root = str(Path(__file__).resolve().parent.parent) + '/'
sys.path.append(project_root)

SPECIES_SYNONYMS = {
    "Poplar": ["Poplar", "Pptrich", "Ptrichocarpa", "ptrich", "ptr"],
    "Sorghum": ["Sorghum", "Sbicolor", "Sbi", "sbi"],
}
SPECIES_ORTHOLOG_FILE = {
    "Poplar": "Ath-Ptr-Orthologs.tsv",
    "Sorghum": "Ath-Sbi-Orthologs.tsv",
}
CONTROL_COLOR = "#8c7999"
TREATMENT_COLOR = "#71798c"
PLASTID_COLOR = "#4e9e6e"
NON_PLASTID_COLOR = "#9b5fa5"
ALL_GENES_COLOR = "#999999"


def find_tmm_file(tmm_dir, species_name):
    synonyms = SPECIES_SYNONYMS.get(species_name, [species_name])
    for file_name in os.listdir(tmm_dir):
        if any(syn.lower() in file_name.lower() for syn in synonyms):
            return os.path.join(tmm_dir, file_name)
    raise FileNotFoundError(f"No TMM file found for {species_name} in {tmm_dir}")


def load_plastid_gene_set(plastid_gene_file):
    with open(plastid_gene_file, 'r') as f:
        return {line.strip().upper() for line in f if line.strip()}


def load_ortholog_map(ortholog_file):
    ortho_df = pa.read_csv(ortholog_file, sep='\t', header=None, usecols=[0, 1])
    ortho_df.columns = ['Ath_ID', 'Species_ID']
    ortho_df = ortho_df.drop_duplicates(subset=['Species_ID'])
    return dict(zip(ortho_df['Species_ID'], ortho_df['Ath_ID']))


def load_classified_tmm(species_name, tmm_dir, orthologs_dir, plastid_gene_file,
                         tissue='Leaf', treatments=('Control', 'FeLim'), exclude_days=('0h', '1h'),
                         include_days=None):
    tmm_file = find_tmm_file(tmm_dir, species_name)
    # CLAUDE 2026-08-12: the TMM tables shipped in projects/*/rnaseq-data are
    # tab-separated and carry a single `condition` column of the form
    # Tissue_Treatment_Timepoint (e.g. Leaf_FeLim_7d) rather than the three
    # separate columns this loader was written against. Sniff the delimiter and
    # split `condition` when present, so the same code reads either layout.
    with open(tmm_file, 'rb') as _fh:
        _head = _fh.read(4096)
    if tmm_file.endswith(('.xz', '.gz', '.bz2')):
        _sep = '\t'
    else:
        _sep = '\t' if _head.split(b'\n')[0].count(b'\t') else ','
    tmm_df = pa.read_csv(tmm_file, sep=_sep)
    if 'condition' in tmm_df.columns and 'tissue' not in tmm_df.columns:
        parts = tmm_df['condition'].str.split('_', n=2, expand=True)
        tmm_df['tissue'], tmm_df['treatment'], tmm_df['time_stamp'] = (
            parts[0], parts[1], parts[2])

    tmm_df = tmm_df[(tmm_df['tissue'] == tissue) & (tmm_df['treatment'].isin(treatments))]
    if exclude_days:
        tmm_df = tmm_df[~tmm_df['time_stamp'].isin(exclude_days)]
    if include_days:
        tmm_df = tmm_df[tmm_df['time_stamp'].isin(include_days)]
    tmm_df = tmm_df[tmm_df['value'] > 0]

    at_plastid_genes = load_plastid_gene_set(plastid_gene_file)
    ortholog_file = os.path.join(orthologs_dir, SPECIES_ORTHOLOG_FILE[species_name])
    ortho_map = load_ortholog_map(ortholog_file)

    def is_plastidial(gene_id):
        at_id = ortho_map.get(gene_id, "")
        at_id_clean = str(at_id).split('.')[0].upper()
        return at_id_clean in at_plastid_genes

    tmm_df = tmm_df.copy()
    tmm_df['is_plastid'] = tmm_df['Gene_ID'].apply(is_plastidial)
    tmm_df['species'] = species_name
    return tmm_df


# Figure 1: count-based abundance density (where plastid genes' expression levels fall
# relative to the rest of the transcriptome), with a rug plot marking individual plastid genes.
def plot_density_with_rug(species_data, species_list, output_path):
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(len(species_list), 1, figsize=(7, 3.2 * len(species_list)),
                              facecolor='white', sharex=False)
    if len(species_list) == 1:
        axes = [axes]

    for ax, species_name in zip(axes, species_list):
        df = species_data[species_name]
        sns.kdeplot(data=df, x='value', hue='treatment', log_scale=True, common_norm=False,
                    palette={'Control': CONTROL_COLOR, 'FeLim': TREATMENT_COLOR},
                    linewidth=2, ax=ax)

        plastid_df = df[df['is_plastid']]
        sns.rugplot(data=plastid_df, x='value', hue='treatment',
                    palette={'Control': CONTROL_COLOR, 'FeLim': TREATMENT_COLOR},
                    height=0.07, alpha=0.35, legend=False, ax=ax)

        ax.set_title(species_name, fontsize=13)
        ax.set_xlabel('TMM-normalized abundance (log scale)', fontsize=11)
        ax.set_ylabel('Density (gene count)', fontsize=11)
        ax.set_facecolor('white')
        if ax.legend_ is not None:
            ax.legend_.set_title(None)

    fig.suptitle('')
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, facecolor='white')
    plt.close(fig)


# Figure 2: mass-weighted abundance density (where most of the *transcript abundance*,
# not just gene count, is concentrated) for plastid genes vs. the full transcriptome.
def plot_mass_weighted_density(species_data, species_list, output_path):
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(len(species_list), 1, figsize=(7, 3.2 * len(species_list)),
                              facecolor='white', sharex=True, sharey=True)
    if len(species_list) == 1:
        axes = [axes]

    for ax, species_name in zip(axes, species_list):
        df = species_data[species_name]

        for treatment, dash in [('Control', False), ('FeLim', True)]:
            all_df = df[df['treatment'] == treatment]

            plastid_df = all_df[all_df['is_plastid']]
            sns.kdeplot(data=plastid_df, x='value', weights='value', log_scale=True,
                        color=PLASTID_COLOR, linestyle='--' if dash else '-',
                        linewidth=2, fill=not dash, alpha=0.25 if not dash else 0.9, ax=ax,
                        label=f'Plastid genes ({treatment})')

            nonplastid_df = all_df[~all_df['is_plastid']]
            sns.kdeplot(data=nonplastid_df, x='value', weights='value', log_scale=True,
                        color=NON_PLASTID_COLOR, linestyle='--' if dash else '-',
                        linewidth=2, fill=not dash, alpha=0.2 if not dash else 0.8, ax=ax,
                        label=f'Non-plastid genes ({treatment})')

        ax.set_title(species_name, fontsize=13)
        ax.set_xlabel('TMM-normalized abundance (log scale)', fontsize=11)
        ax.set_ylabel('Density (mass-weighted)', fontsize=11)
        ax.set_facecolor('white')
        ax.set_ylim(0, 0.7)
        ax.legend(fontsize=8, loc='upper left')

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, facecolor='white')
    plt.close(fig)


# Figure 3: per-day count-based KDE grid — species as rows, days as columns.
# Each panel shows plastid vs non-plastid gene density (log-scale x), with
# solid = Control and dashed = FeLim.
def plot_density_per_day(species_data, species_list, output_path,
                          day_order=('2d', '4d', '7d', '14d', '21d')):
    days = list(day_order)
    n_rows, n_cols = len(species_list), len(days)

    sns.set_style('whitegrid')
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.8 * n_cols, 2.6 * n_rows),
        facecolor='white',
        sharex=True, sharey=True,
    )
    if n_rows == 1:
        axes = [axes]

    for row_idx, species_name in enumerate(species_list):
        df = species_data[species_name]
        for col_idx, day in enumerate(days):
            ax = axes[row_idx][col_idx]
            day_df = df[df['time_stamp'] == day]

            for treatment, ls in [('Control', '-'), ('FeLim', '--')]:
                trmt_df = day_df[day_df['treatment'] == treatment]
                if trmt_df.empty:
                    continue

                plastid_df    = trmt_df[trmt_df['is_plastid']]
                nonplastid_df = trmt_df[~trmt_df['is_plastid']]

                label_p  = f'Plastid ({treatment})'    if row_idx == 0 and col_idx == 0 else None
                label_np = f'Non-plastid ({treatment})' if row_idx == 0 and col_idx == 0 else None

                if len(plastid_df) > 5:
                    sns.kdeplot(data=plastid_df, x='value', log_scale=True,
                                color=PLASTID_COLOR, linestyle=ls, linewidth=1.6,
                                ax=ax, label=label_p)
                if len(nonplastid_df) > 5:
                    sns.kdeplot(data=nonplastid_df, x='value', log_scale=True,
                                color=NON_PLASTID_COLOR, linestyle=ls, linewidth=1.6,
                                ax=ax, label=label_np)

            ax.set_facecolor('white')
            ax.set_ylim(0, 0.7)
            ax.tick_params(labelsize=8)
            if col_idx == 0:
                ax.set_ylabel(species_name, fontsize=10, fontweight='bold')
            else:
                ax.set_ylabel('')
            if row_idx == 0:
                ax.set_title(day, fontsize=10, fontweight='bold')
            if row_idx == n_rows - 1:
                ax.set_xlabel('TMM abundance', fontsize=8)
            else:
                ax.set_xlabel('')

    # Single shared legend from the first panel's handles
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', ncol=4,
               fontsize=8, frameon=True, bbox_to_anchor=(0.98, 0.56),
               labelspacing=0.2, handletextpad=0.4, handlelength=1.5, borderpad=0.3)

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(output_path, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate abundance distribution figures highlighting plastid-gene location/contribution."
    )
    parser.add_argument("--species", nargs='+', default=["Poplar", "Sorghum"])
    parser.add_argument("--tmm-dir", default=os.path.join(project_root, "data", "RNAseq", "tmm"))
    parser.add_argument("--orthologs-dir", default=os.path.join(project_root, "data", "orthologs"))
    parser.add_argument("--plastid-gene-file",
                         default=os.path.join(project_root, "data", "plastid_proteome",
                                               "compiled_arabidopsis_plastid_genes.txt"))
    parser.add_argument("--tissue", default="Leaf")
    parser.add_argument("--treatments", nargs='+', default=["Control", "FeLim"])
    parser.add_argument("--exclude-days", nargs='*', default=["0h", "1h"])
    parser.add_argument("--include-days", nargs='*', default=None,
                         help="If set, only use data from these time points (e.g. --include-days 7d)")
    parser.add_argument("--output-density-rug", default="plastid_abundance_density_rug.png")
    parser.add_argument("--output-mass-weighted", default="plastid_abundance_mass_weighted.png")
    parser.add_argument("--output-mass-weighted-pdf", default="plastid_abundance_mass_weighted.pdf")
    parser.add_argument("--output-per-day", default="plastid_abundance_per_day.png")
    parser.add_argument("--output-per-day-pdf", default="plastid_abundance_per_day.pdf")
    args = parser.parse_args()

    species_data = {}
    for species_name in args.species:
        species_data[species_name] = load_classified_tmm(
            species_name, args.tmm_dir, args.orthologs_dir, args.plastid_gene_file,
            tissue=args.tissue, treatments=tuple(args.treatments),
            exclude_days=tuple(args.exclude_days) if args.exclude_days else (),
            include_days=tuple(args.include_days) if args.include_days else None)

    plot_density_with_rug(species_data, args.species, args.output_density_rug)
    print(f"Saved {args.output_density_rug}")

    plot_mass_weighted_density(species_data, args.species, args.output_mass_weighted)
    print(f"Saved {args.output_mass_weighted}")
    plot_mass_weighted_density(species_data, args.species, args.output_mass_weighted_pdf)
    print(f"Saved {args.output_mass_weighted_pdf}")

    plot_density_per_day(species_data, args.species, args.output_per_day)
    print(f"Saved {args.output_per_day}")
    plot_density_per_day(species_data, args.species, args.output_per_day_pdf)
    print(f"Saved {args.output_per_day_pdf}")
