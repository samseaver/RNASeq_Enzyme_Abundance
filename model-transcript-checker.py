#!/usr/bin/env python
import cobra
import csv
import pandas as pd
import re

from cobrakbase.core.kbase_object_factory import KBaseObjectFactory

# --- USER INPUTS ---
model_path = "projects/brave/inputs/Sbicolor-v5.1-reconstruction.json"
transcript_path = "projects/brave/rnaseq-data/tmm/average_genecounts_long_Sbicolor.tsv"
transcript_id_col = "Geneid"
# -------------------

def check_id_mapping():
    print(f"Loading model from: {model_path}")
    try:
        # Auto-detect model format (SBML, JSON, YAML)
        if model_path.endswith('.json'):
            KBOF = KBaseObjectFactory()
            model = KBOF.build_object_from_file(model_path, "KBaseFBA.FBAModel")
        elif model_path.endswith('.yaml'):
            model = cobra.io.load_yaml_model(model_path)
        else:
            model = cobra.io.read_sbml_model(model_path)
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    print(f"Loading transcripts from: {transcript_path}")
    try:
        # Read the first few bytes to guess the separator
        with open(transcript_path, 'r') as f:
            # 'sniff' the first 1024 bytes to find the dialect
            dialect = csv.Sniffer().sniff(f.read(1024))
    
        # Now read it with the detected delimiter using the fast C engine (default)
        df = pd.read_csv(transcript_path, sep=dialect.delimiter)
    except Exception as e:
        print(f"Error loading transcript data: {e}")
        return

    # 1. Get Sets of IDs
    model_gene_ids = {g.id for g in model.genes}
    transcript_ids = set(df[transcript_id_col].astype(str))

    print(f"\n--- Initial State ---")
    print(f"Total Model Genes: {len(model_gene_ids)}")
    print(f"Total Transcript IDs: {len(transcript_ids)}")
    
    # Check raw overlap
    raw_overlap = model_gene_ids.intersection(transcript_ids)
    print(f"Direct Match Count: {len(raw_overlap)}")

    # 2. Simulate the Clean-up
    # Regex to remove suffix: matches a dot, one or more digits, dot, 'p' at the end
    regex_pattern = r'\.\d+\.p$' 
    
    cleaned_model_map = {}
    for original_id in model_gene_ids:
        clean_id = re.sub(regex_pattern, '', original_id)
        cleaned_model_map[original_id] = clean_id

    cleaned_model_ids = set(cleaned_model_map.values())
    
    # Check overlap after cleaning
    cleaned_overlap = cleaned_model_ids.intersection(transcript_ids)
    
    print(f"\n--- After Applying Regex (r'{regex_pattern}') ---")
    print(f"Matches found: {len(cleaned_overlap)}")
    
    if len(cleaned_overlap) > 0:
        percent_match = (len(cleaned_overlap) / len(cleaned_model_ids)) * 100
        print(f"Percentage of model genes covered by transcript data: {percent_match:.2f}%")
        
    else:
        print("\nWARNING: Still found 0 matches after cleaning. Check regex or ID formats.")

    # 3. Mismatch Analysis (New Section)
    unmatched_originals = []
    
    for original, clean in cleaned_model_map.items():
        if clean not in transcript_ids:
            unmatched_originals.append(original)

    print(f"\n--- Mismatch Analysis ---")
    print(f"Unmatched Genes: {len(unmatched_originals)}")

    if len(unmatched_originals) > 0:
        print("\n--- Sample of UNMATCHED Model IDs (Original -> Cleaned) ---")
        # Print first 20 to spot patterns
        for x in unmatched_originals[:20]:
            clean_version = cleaned_model_map[x]
            print(f"Original: {x:<25} -> Cleaned: {clean_version}")

        print("\n--- Diagnosis Hint ---")
        print("1. If 'Cleaned' looks correct (e.g., Sobic.001G...), then these genes are missing from your transcript file.")
        print("2. If 'Cleaned' looks wrong (still has suffixes), the regex needs adjustment.")
        
        # Check for partial matches
        print("\n--- Checking for partial matches in transcript data ---")
        if unmatched_originals:
            sample_fail = cleaned_model_map[unmatched_originals[0]]
            print(f"Searching transcript IDs for anything resembling '{sample_fail}'...")
            
            transcript_list = list(transcript_ids)
            similar = [t for t in transcript_list[:5000] if sample_fail in t] # quick scan of first 5000
            
            if similar:
                print(f"Found similar transcript IDs: {similar[:5]}")
            else:
                print(f"No similar IDs found in the first 5000 transcript records for {sample_fail}.")

if __name__ == "__main__":
    check_id_mapping()