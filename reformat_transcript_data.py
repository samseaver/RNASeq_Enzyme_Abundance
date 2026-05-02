#!/usr/bin/env python
import pandas as pd
import sys

# --- CONFIGURATION ---
INPUT_FILE = 'Poplar_raw_genes_tmm_mean_std.csv'  # Replace with your actual file name
OUTPUT_FILE = 'Poplar_raw_genes_tmm_mean.tsv'

# 1. Load the data
# Assumes the input is comma-separated based on your snippet
try:
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} rows from {INPUT_FILE}")
except FileNotFoundError:
    print(f"Error: Could not find {INPUT_FILE}")
    sys.exit(1)

# 2. Create the new concatenated 'condition' column
# We join them with underscores: Leaf_Control_0h
df['condition'] = df['tissue'] + '_' + df['treatment'] + '_' + df['time_stamp']

# 3. Select and Reorder columns
# Dropping the original 3 columns and keeping the rest
cols_to_keep = ['Gene_ID', 'condition', 'value', 'value_std', 'value_log', 'value_std_log']
df_out = df[cols_to_keep]

# 4. Save as Tab-Separated Values (TSV)
df_out.to_csv(OUTPUT_FILE, sep='\t', index=False)
print(f"Saved formatted data to {OUTPUT_FILE}")

# Optional: Print the first few lines to verify
print("\nSample Output:")
print(df_out.head())