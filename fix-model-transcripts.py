#!/usr/bin/env python
import json
import re
import sys

# --- USER INPUTS ---
input_model_path = "Sbicolor-v5.1-reconstruction.json"
output_model_path = "projects/brave/inputs/Sbicolor-v5.1-reconstruction_fixed.json"

# Sample Regex to remove suffix: matches a dot, one or more digits, dot, 'p' at the end
TRANSCRIPT_REGEX = r'\.\d+\.p$'
# -------------------

def fix_model_json():
    print(f"Loading raw JSON from: {input_model_path}")
    
    try:
        with open(input_model_path, 'r') as f:
            model_data = json.load(f)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return
    
    # Track changes
    gene_map = {} # old_id_suffix -> new_id_clean

    # KBase models link reactions to genes via 'modelReactionProteins' -> 'feature_refs'
    refs_updated = 0
            
    for rxn in model_data['modelreactions']:
        # Check if this reaction has protein data
        if 'modelReactionProteins' in rxn:
            for protein in rxn['modelReactionProteins']:
                # Check for subunits
                if 'modelReactionProteinSubunits' in protein:
                    for subunit in protein['modelReactionProteinSubunits']:
                        # 'feature_refs' is the list of strings pointing to genes
                        if 'feature_refs' in subunit:
                            new_refs = []
                            for old_ref in subunit['feature_refs']:
                                # Ref format is usually: "~/genomes/GUID/features/GeneID"
                                # We need to check if the END of the ref matches our old ID
                                    
                                matched = False
                                # Apply regex
                                new_ref = re.sub(TRANSCRIPT_REGEX, '', old_ref)
                                if(new_ref != old_ref):
                                    gene_map[old_ref] = new_ref
                                    new_refs.append(new_ref)
                                    refs_updated += 1
                                    matched = True

                                if not matched:
                                    # If no change needed, keep original
                                    new_refs.append(old_ref)
                                
                            # Assign the updated list back
                            subunit['feature_refs'] = new_refs
    
    print(f"Updated {refs_updated} feature references.")

    # ---------------------------------------------------------
    # STEP 3: Save
    # ---------------------------------------------------------
    print(f"Saving to: {output_model_path}")
    try:
        with open(output_model_path, 'w') as f:
            json.dump(model_data, f, indent=4)
        print("Success! Model saved.")
    except Exception as e:
        print(f"Error saving JSON: {e}")

if __name__ == "__main__":
    fix_model_json()