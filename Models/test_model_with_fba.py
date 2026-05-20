#!/usr/bin/env python3
import argparse
import sys
import json
import cobra
from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
from cobra.flux_analysis import flux_variability_analysis as fva

def load_model_auto(model_path):
    """Auto-detects model format and loads it into a cobra model object."""
    if model_path.endswith('.json'):
        # Peek inside to distinguish between KBase JSON and standard COBRA JSON
        try:
            with open(model_path, 'r') as f:
                data = json.load(f)
            
            if 'modelreactions' in data:
                print(" -> Detected KBase JSON format. Routing to KBaseObjectFactory...")
                KBOF = KBaseObjectFactory()
                return KBOF.build_object_from_file(model_path, "KBaseFBA.FBAModel")
            else:
                print(" -> Detected Standard COBRA JSON format...")
                return cobra.io.load_json_model(model_path)
        except Exception as e:
            sys.exit(f"Failed to parse JSON file: {e}")
            
    elif model_path.endswith('.yaml') or model_path.endswith('.yml'):
        print(" -> Detected YAML format...")
        return cobra.io.load_yaml_model(model_path)
    else:
        print(" -> Detected SBML/XML format...")
        return cobra.io.read_sbml_model(model_path)

def test_model(model_path, objective='bio1', run_fva=False):
    print(f"Loading model from: {model_path}...")
    try:
        model = load_model_auto(model_path)
        print(f" -> Successfully loaded model with {len(model.reactions)} reactions and {len(model.metabolites)} metabolites.\n")
    except Exception as e:
        sys.exit(f"Error loading model: {e}")
    
    print(f"Optimizing model for objective: {objective}...")
    try:
        solution = model.optimize(objective)
        print(f" -> Optimization Status: {solution.status}")
        print(f" -> Objective Value (Growth): {solution.objective_value}\n")
    except Exception as e:
        print(f"Optimization failed: {e}")
        return
        
    if run_fva:
        print("Running Flux Variability Analysis (FVA)...")
        print("(Note: This may take a while depending on the network size.)")
        try:
            fva_result = fva(model, fraction_of_optimum=1.0, processes=1)
            print("\n=== FVA Results ===")
            print(fva_result)
        except Exception as e:
            print(f"FVA failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test loading and optimizing a metabolic model (Auto-detects format).")
    parser.add_argument("-m", "--model", required=True, help="Path to the input model file (.json, .yaml, .xml).")
    parser.add_argument("-o", "--objective", default="bio1", help="Objective function to optimize (Default: bio1).")
    parser.add_argument("--fva", action="store_true", help="Include this flag to run Flux Variability Analysis (FVA).")
    
    args = parser.parse_args()
    
    test_model(args.model, args.objective, args.fva)