#!/usr/bin/env python3
import pandas as pd
import os
import sys

# --- HARD-CODED CONFIGURATION ---
at_chloro_file = 'Full_AT_CHLORO_2019.tsv'
ppdb_file = 'PPDB_Plastid_Proteome.tsv'
output_file = 'compiled_arabidopsis_plastid_genes.txt'

def clean_at_id(gene_id):
    """Strips transcript splice-variant suffixes (e.g., .1, .2) to return base AT gene ID."""
    return str(gene_id).split('.')[0].strip().upper()

def main():
    print("=== Internal Pipeline: Building Arabidopsis Plastid Database ===")
    
    # 1. Process AT_CHLORO Database
    if not os.path.exists(at_chloro_file):
        sys.exit(f"Error: Could not find AT_CHLORO file at {at_chloro_file}")
        
    print(f" -> Loading and filtering AT_CHLORO...")
    at_chloro_df = pd.read_csv(at_chloro_file, sep='\t')
    
    # Collapse all columns to check for sub-compartment tags (THY, STR, ENV)
    row_strings = at_chloro_df.astype(str).apply(' '.join, axis=1)
    mask = row_strings.str.contains('THY|STR|ENV', case=True, regex=True)
    filtered_at_chloro_df = at_chloro_df[mask]
    
    # FIX: Based on your grep, the Gene ID is in the SECOND column (index 1)
    # Column 0 appears to be a numeric index
    at_chloro_genes = set(filtered_at_chloro_df.iloc[:, 1].dropna().apply(clean_at_id))
    print(f"    - Extracted {len(at_chloro_genes)} base genes from filtered AT_CHLORO.")

    # 2. Process PPDB Database
    if not os.path.exists(ppdb_file):
        sys.exit(f"Error: Could not find PPDB file at {ppdb_file}")

    print(f" -> Loading PPDB...")
    # Assuming PPDB starts with the ID in the first column
    ppdb_df = pd.read_csv(ppdb_file, sep='\t')
    ppdb_genes = set(ppdb_df.iloc[:, 0].dropna().apply(clean_at_id))
    print(f"    - Extracted {len(ppdb_genes)} base genes from PPDB.")

    # 3. Merge and Write
    at_plastid_union = at_chloro_genes.union(ppdb_genes)
    
    # Final filter: Ensure we only keep IDs starting with 'AT' 
    # This prevents any stray index numbers from leaking into the final list
    final_genes = {g for g in at_plastid_union if g.startswith('AT')}
    
    print(f"\n -> Total unique Plastid genes (AT-prefixed only): {len(final_genes)}")
    
    try:
        with open(output_file, 'w') as out:
            for gene in sorted(final_genes):
                out.write(f"{gene}\n")
        print(f"SUCCESS: Master list saved to {output_file}")
    except Exception as e:
        print(f"CRITICAL ERROR: Could not save output file: {e}")

if __name__ == "__main__":
    main()
