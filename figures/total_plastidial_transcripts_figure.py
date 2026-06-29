#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# --- Configuration ---
input_file = "../projects/qpsi-plastidial/integration_results/Poplar_plastid_transcript_totals.tsv"
output_plot = "Poplar_plastid_transcripts_plot.png"

if not os.path.exists(input_file):
    print(f"Error: Could not find {input_file}")
    sys.exit(1)

print(f"Reading data from {input_file}...")
df = pd.read_csv(input_file, sep='\t')

# Extract treatment group (e.g., 'Leaf_Control') and time point (e.g., '14d')
df['Treatment'] = df['condition'].apply(lambda x: '_'.join(x.split('_')[:-1]))
df['Time_Label'] = df['condition'].apply(lambda x: x.split('_')[-1])

# Convert time labels to numeric hours for proper chronological sorting and plotting
def time_to_hours(t):
    if t.endswith('h'):
        return float(t.replace('h', ''))
    elif t.endswith('d'):
        return float(t.replace('d', '')) * 24
    return 0

df['Hours'] = df['Time_Label'].apply(time_to_hours)

# Sort the dataframe to ensure the lines connect chronologically
df = df.sort_values(by=['Treatment', 'Hours'])

# --- Plotting ---
print("Generating plot...")
plt.figure(figsize=(10, 6))

for treatment in df['Treatment'].unique():
    subset = df[df['Treatment'] == treatment]
    # Plot against Hours to ensure the x-axis scale accurately represents the time gaps
    plt.plot(subset['Hours'], subset['total_plastid_transcripts'], marker='o', linewidth=2, label=treatment)

# Map the hours back to your original labels (e.g., '0h', '2d') for a cleaner x-axis
all_hours = sorted(df['Hours'].unique())
hour_to_label = dict(zip(df['Hours'], df['Time_Label']))
labels = [hour_to_label[h] for h in all_hours]

plt.xticks(all_hours, labels)
plt.xlabel('Time')
plt.ylabel('Total Plastid Transcripts (Sum of Abundance)')
plt.title('Change in Total Plastidial Transcripts Over Time')
plt.legend(title='Treatment Group')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

plt.savefig(output_plot, dpi=300)
print(f"Plot saved successfully to {output_plot}")