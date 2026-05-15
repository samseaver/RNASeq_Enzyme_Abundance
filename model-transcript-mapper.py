#!/usr/bin/env python3
import cobra
import csv
import pandas as pd
import re
import argparse
import sys
import json
import os
from collections import Counter

def load_transcripts(transcript_path, transcript_id_col):
    """Loads the transcript file and returns a set of unique Gene IDs."""
    print(f"Loading transcripts from: {transcript_path}")
    try:
        with open(transcript_path, 'r') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
        df = pd.read_csv(transcript_path, sep=dialect.delimiter)
        
        if transcript_id_col not in df.columns:
            print(f"Error: Column '{transcript_id_col}' not found in transcript file.")
            print(f"Available columns: {list(df.columns)}")
            sys.exit(1)
            
        transcript_ids = set(df[transcript_id_col].astype(str).str.strip())
        print(f" -> Found {len(transcript_ids)} unique transcript IDs.\n")
        return transcript_ids
    except Exception as e:
        print(f"Error reading transcript file: {e}")
        sys.exit(1)

def extract_model_genes(model_path):
    """Extracts raw gene IDs depending on model format."""
    print(f"Loading model from: {model_path}")
    model_data = None
    model_genes = []
    model_type = None

    try:
        if model_path.endswith('.json'):
            model_type = 'kbase_json'
            with open(model_path, 'r') as f:
                model_data = json.load(f)
            # KBase JSON extraction
            for rxn in model_data.get('modelreactions', []):
                for prot in rxn.get('modelReactionProteins', []):
                    for sub in prot.get('modelReactionProteinSubunits', []):
                        for ref in sub.get('feature_refs', []):
                            gene_id = ref.split('/')[-1]
                            if gene_id not in model_genes:
                                model_genes.append(gene_id)
        else:
            model_type = 'cobra'
            if model_path.endswith('.yaml') or model_path.endswith('.yml'):
                model_data = cobra.io.load_yaml_model(model_path)
            else:
                model_data = cobra.io.read_sbml_model(model_path)
            model_genes = [g.id for g in model_data.genes]
            
        print(f" -> Found {len(model_genes)} unique genes in the metabolic model.\n")
        return model_type, model_data, model_genes
    except Exception as e:
        print(f"Error loading model: {e}")
        sys.exit(1)

def predict_regex(unmatched_model_genes, transcript_ids):
    """Attempts to guess the regex required by finding common prefixes."""
    print("\n=== AUTO-DETECTING REGEX FIXES ===")
    
    sample_model = list(unmatched_model_genes)[:5]
    sample_transcript = list(transcript_ids)[:5]
    
    # Updated Labels
    print(f"Sample Unmatched Model IDs:  {sample_model}")
    print(f"Sample Valid Transcript IDs: {sample_transcript}")
    
    # New Explicit Warning
    print("\n[INFO] Note: If the unmatched model IDs look structurally completely different")
    print("from your valid transcript IDs (e.g., organelle genes vs nuclear genes), those")
    print("missing transcripts may simply not be present in your RNA-Seq dataset at all.")
    
    predicted_patterns = Counter()
    
    # Check the first 500 unmatched genes to guess the pattern (for speed)
    for m_id in list(unmatched_model_genes)[:500]:
        for t_id in transcript_ids:
            # If the model ID starts with the transcript ID, the rest is the "junk" suffix
            if m_id.startswith(t_id) and len(m_id) > len(t_id):
                raw_suffix = m_id[len(t_id):]
                
                # Convert the raw suffix into a generalized regex
                # Example: '.CDS.1' -> '\.CDS\.\d+'
                generalized = raw_suffix.replace('.', r'\.')
                generalized = re.sub(r'\d+', r'\\d+', generalized)
                pattern = generalized + r'$'
                
                predicted_patterns[pattern] += 1
                break # Move to next model gene once a root is found
                
    if predicted_patterns:
        print("\n[!] SUCCESS: Possible suffix patterns detected!")
        print("Try running the script again with one of these suggested arguments:")
        for pattern, count in predicted_patterns.most_common(3):
            print(f"    -r '{pattern}'    (Successfully mapped {count} sample IDs)")
    else:
        print("\n[?] Could not auto-detect a simple suffix difference.")
        print("The IDs might be completely different formats (e.g., UniProt vs AGI),")
        print("or they are simply the missing transcripts mentioned in the INFO block above.")
    print("==================================\n")

def main():
    parser = argparse.ArgumentParser(description="Check and Fix mappings between model genes and transcript IDs.")
    parser.add_argument("-m", "--model", required=True, help="Path to the model file (.json, .yaml, or .xml)")
    parser.add_argument("-t", "--transcript", required=True, help="Path to the transcript file (TSV/CSV)")
    parser.add_argument("-c", "--column", default="Gene_ID", help="Name of column containing transcript IDs (Default: Gene_ID)")
    parser.add_argument("-r", "--regex", nargs='*', default=[], help="One or more regex patterns to strip from model IDs. Evaluated in order.")
    parser.add_argument("-o", "--output", help="Optional. If provided, applies the successful regex fixes and saves a new model.")
    
    args = parser.parse_args()

    # 1. Load Data
    transcript_ids = load_transcripts(args.transcript, args.column)
    model_type, model_data, model_genes = extract_model_genes(args.model)

    # 2. Map and Check Mismatches
    gene_map = {} 
    
    direct_matches = 0
    regex_matches = {pattern: 0 for pattern in args.regex}
    unmatched_originals = []
    
    for original in model_genes:
        # Check Direct Match
        if original in transcript_ids:
            gene_map[original] = original
            direct_matches += 1
            continue
            
        # Check Regex Transformations
        matched = False
        for pattern in args.regex:
            clean_guess = re.sub(pattern, '', original)
            if clean_guess != original and clean_guess in transcript_ids:
                gene_map[original] = clean_guess
                regex_matches[pattern] += 1
                matched = True
                break
                
        if not matched:
            unmatched_originals.append(original)

    # 3. Print Report
    total_matched = direct_matches + sum(regex_matches.values())
    
    print("=== MAPPING REPORT ===")
    print(f"Total Model Genes:    {len(model_genes)}")
    print(f"Direct Matches:       {direct_matches}")
    
    for pattern in args.regex:
        print(f"Matches via Regex '{pattern}': {regex_matches[pattern]}")
        
    print(f"Total Matched:        {total_matched}")
    print(f"Still Unmatched:      {len(unmatched_originals)}\n")

    # 4. Auto-Predict Regex if there are still unmatched genes
    if len(unmatched_originals) > 0:
        predict_regex(unmatched_originals, transcript_ids)

    # 5. Save Fixed Model (If requested)
    if args.output:
        if total_matched == 0:
            print("Aborting save: No genes matched. Check your regex or transcript file.")
            sys.exit(1)
            
        print(f"=== APPLYING FIXES AND SAVING ===")
        print(f"Saving fixed model to: {args.output}")
        
        if model_type == 'kbase_json':
            refs_updated = 0
            for rxn in model_data['modelreactions']:
                if 'modelReactionProteins' in rxn:
                    for protein in rxn['modelReactionProteins']:
                        if 'modelReactionProteinSubunits' in protein:
                            for subunit in protein['modelReactionProteinSubunits']:
                                new_refs = []
                                for ref in subunit.get('feature_refs', []):
                                    base_ref = ref.rsplit('/', 1)[0]
                                    old_id = ref.split('/')[-1]
                                    
                                    if old_id in gene_map and gene_map[old_id] != old_id:
                                        new_id = gene_map[old_id]
                                        new_refs.append(f"{base_ref}/{new_id}")
                                        refs_updated += 1
                                    else:
                                        new_refs.append(ref)
                                subunit['feature_refs'] = new_refs
                                
            with open(args.output, 'w') as f:
                json.dump(model_data, f, indent=4)
            print(f"Success! Updated {refs_updated} KBase feature references.")
            
        elif model_type == 'cobra':
            genes_updated = 0
            for gene in model_data.genes:
                if gene.id in gene_map and gene_map[gene.id] != gene.id:
                    gene.id = gene_map[gene.id]
                    genes_updated += 1
                    
            if args.output.endswith('.yaml') or args.output.endswith('.yml'):
                cobra.io.save_yaml_model(model_data, args.output)
            else:
                cobra.io.write_sbml_model(model_data, args.output)
            print(f"Success! Updated {genes_updated} Cobra gene IDs.")

if __name__ == "__main__":
    main()