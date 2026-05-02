#!/usr/bin/env python
import pandas as pd
import os

# --- Configuration: File Paths ---
project_root = "." # Update this if running from a different directory
species,short_species = "Sorghum","Sbi"
species,short_species = "Poplar","Ptr"

ortholog_file = os.path.join(project_root, f"data/orthologs/Ath-{short_species}-Orthologs.tsv")
tmm_file = os.path.join(project_root, f"projects/qpsi-plastidial/rnaseq-data/tmm/{species}_raw_genes_tmm_mean.tsv")
rxn_abundance_file = os.path.join(project_root, f"projects/qpsi-plastidial/integration_results/{species}_objective_abundance.tsv")

at_chloro_file = os.path.join(project_root, 'data/plastid_proteome/Full_AT_CHLORO_2019.tsv')
ppdb_file = os.path.join(project_root, 'data/plastid_proteome/PPDB_Plastid_Proteome.tsv')

# Outputs
fractions_output_file = os.path.join(project_root, f"projects/qpsi-plastidial/integration_results/{species}_rxn_molar_fractions.tsv")
totals_output_file = os.path.join(project_root, f"projects/qpsi-plastidial/integration_results/{species}_plastid_transcript_totals.tsv")

def clean_at_id(gene_id):
    """Strips transcript suffixes (e.g., .1, .2) to return base AT gene ID."""
    return str(gene_id).split('.')[0].strip().upper()

print("1. Building Arabidopsis Plastid Proteome Union...")
# Explicit usecols ensures we don't load extra garbage from hidden columns
at_chloro_df = pd.read_csv(at_chloro_file, sep='\t')
row_strings = at_chloro_df.astype(str).apply(' '.join, axis=1)
mask = row_strings.str.contains('THY|STR|ENV', case=True, regex=True)
filtered_at_chloro_df = at_chloro_df[mask]
at_chloro_genes = set(filtered_at_chloro_df.iloc[:, 0].dropna().apply(clean_at_id))

ppdb_df = pd.read_csv(ppdb_file, sep='\t', usecols=[0]) 
ppdb_genes = set(ppdb_df.iloc[:, 0].dropna().apply(clean_at_id))

at_plastid_union = at_chloro_genes.union(ppdb_genes)
print(f"   -> Found {len(at_plastid_union)} unique Arabidopsis plastid genes.")


print("2. Mapping to Poplar Plastid Proteome...")
# CRITICAL FIX: usecols=[0,1] prevents pandas from turning extra columns into an index
ortho_df = pd.read_csv(ortholog_file, sep='\t', header=None, usecols=[0, 1], names=['AT_gene', 'Ptri_gene'])
ortho_df['AT_base'] = ortho_df['AT_gene'].apply(clean_at_id)

poplar_plastid_genes = set(
    ortho_df[ortho_df['AT_base'].isin(at_plastid_union)]['Ptri_gene'].dropna().apply(lambda x: str(x).strip())
)
print(f"   -> Mapped to {len(poplar_plastid_genes)} unique Poplar plastid genes.")


print("3. Calculating Total Plastid Transcript Abundance per Condition...")
# Read TMM data
tmm_df = pd.read_csv(tmm_file, sep='\t')

# CRITICAL FIX: The file has 6 columns. We just extract the 3 we need by their existing names!
tmm_df = tmm_df[['Gene_ID', 'condition', 'value']]

# Filter TMM data to ONLY include genes in the Poplar plastid proteome
plastid_tmm_df = tmm_df[tmm_df['Gene_ID'].isin(poplar_plastid_genes)]

# Group by condition and sum the values to get the denominator for each timepoint
plastid_totals_by_condition = plastid_tmm_df.groupby('condition')['value'].sum().to_dict()

# Save the totals to a separate file
print(f"   -> Saving condition totals to {totals_output_file}...")
totals_df = pd.DataFrame(list(plastid_totals_by_condition.items()), columns=['condition', 'total_plastid_transcripts'])
totals_df.to_csv(totals_output_file, sep='\t', index=False)

for cond, total in plastid_totals_by_condition.items():
    print(f"      - {cond}: {total:.2f}")


print("4. Calculating Molar Fractions for Reactions...")
rxn_df = pd.read_csv(rxn_abundance_file, sep='\t')

def calculate_fraction(row):
    cond = row['condition']
    score = row['reaction_score']
    total = plastid_totals_by_condition.get(cond, 0)
    
    if total > 0:
        return score / total
    return 0.0

rxn_df['molar_fraction'] = rxn_df.apply(calculate_fraction, axis=1)

print(f"5. Saving results to {fractions_output_file}...")
final_df = rxn_df[['condition', 'rxn_ID', 'reaction_score', 'molar_fraction']]
final_df.to_csv(fractions_output_file, sep='\t', index=False)

print("Done!")