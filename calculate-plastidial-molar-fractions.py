#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import os

def load_plastid_genes(plastid_file):
    """Reads the compiled single-column text file of Arabidopsis plastid genes."""
    print(f"Loading master Arabidopsis plastid list: {plastid_file}")
    try:
        with open(plastid_file, 'r') as f:
            genes = {line.strip().upper() for line in f if line.strip()}
        print(f" -> Found {len(genes)} unique Arabidopsis plastid genes.")
        return genes
    except Exception as e:
        sys.exit(f"Error reading plastid list: {e}")

def calculate_molar_fractions(tmm_file, rxn_file, plastid_file, ortholog_file, 
                              id_col, exp_col, val_col, out_rxn):
    
    # 1. Load the Master Arabidopsis Plastid List
    at_plastid_genes = load_plastid_genes(plastid_file)

    # 2. Load TMM Data
    print(f"\nLoading expression data: {tmm_file}")
    tmm_df = pd.read_csv(tmm_file, sep='\t')

    # 3. Handle Ortholog Mapping (Ath Column 0, Species Column 1)
    if ortholog_file and os.path.exists(ortholog_file):
        print(f"Loading ortholog mapping: {ortholog_file}")
        # header=None based on your head snippet; usecols 0 and 1
        ortho_df = pd.read_csv(ortholog_file, sep='\t', header=None, usecols=[0, 1])
        ortho_df.columns = ['Ath_ID', 'Species_ID']
        
        # Filter out duplicates to ensure 1:1 or N:1 mapping
        ortho_df = ortho_df.drop_duplicates(subset=['Species_ID'])
        ortho_map = dict(zip(ortho_df['Species_ID'], ortho_df['Ath_ID']))
        
        def is_plastidial(gene_id):
            at_id = ortho_map.get(gene_id, "")
            at_id_clean = str(at_id).split('.')[0].upper()
            return at_id_clean in at_plastid_genes
            
        print(f" -> Filtering TMM data via Ath-orthology...")
        plastid_tmm_df = tmm_df[tmm_df[id_col].apply(is_plastidial)].copy()
    else:
        print(" -> No ortholog file. Filtering TMM data directly against Ath IDs...")
        plastid_tmm_df = tmm_df[tmm_df[id_col].str.split('.').str[0].str.upper().isin(at_plastid_genes)].copy()

    if plastid_tmm_df.empty:
        sys.exit("Error: No genes in TMM data matched the plastid database. Check ID formats.")

    # 4. Calculate and Save Totals
    print("\nCalculating total plastidial transcripts per condition...")
    plastid_totals = plastid_tmm_df.groupby(exp_col)[val_col].sum().to_dict()
    
    # Auto-generate totals filename based on output name
    totals_file = rxn_file.replace('_reaction_score.tsv', '_plastid_transcript_totals.tsv')
    print(f" -> Saving experiment totals to: {totals_file}")
    
    totals_list = []
    for cond, total in sorted(plastid_totals.items()):
        # print(f"      - {cond}: {total:.2f}")
        totals_list.append({exp_col: cond, 'total_plastid_transcripts': total})
    
    pd.DataFrame(totals_list).to_csv(totals_file, sep='\t', index=False)

    # 5. Calculate Relative Reaction Scores
    print(f"\nLoading reaction abundances: {rxn_file}")
    rxn_df = pd.read_csv(rxn_file, sep='\t')
    
    def get_fraction(row):
        total = plastid_totals.get(row[exp_col], 0)
        return row['reaction_score'] / total if total > 0 else 0.0

    rxn_df['relative_reaction_score'] = rxn_df.apply(get_fraction, axis=1)

    # 6. Save Final Reaction Scores
    rxn_df.to_csv(out_rxn, sep='\t', index=False)
    print(f"Success! Saved results to {out_rxn}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate molar fractions using Arabidopsis orthology.")
    parser.add_argument("-t", "--tmm", required=True, help="TMM expression data (TSV).")
    parser.add_argument("-r", "--rxn", required=True, help="Reaction abundance data (TSV).")
    parser.add_argument("-g", "--genes", default="data/plastid_proteome/compiled_arabidopsis_plastid_genes.txt", help="Master Ath plastid list.")
    parser.add_argument("--orthologs", help="Path to Ath-Species ortholog mapping TSV.")
    parser.add_argument("--id-col", default="Gene_ID", help="Gene ID column in TMM.")
    parser.add_argument("--exp-col", default="condition", help="Condition column.")
    parser.add_argument("--val-col", default="value", help="Expression value column.")
    parser.add_argument("-o", "--output", required=True, help="Output file for relative scores.")
    
    args = parser.parse_args()
    
    calculate_molar_fractions(args.tmm, args.rxn, args.genes, args.orthologs, 
                              args.id_col, args.exp_col, args.val_col, args.output)