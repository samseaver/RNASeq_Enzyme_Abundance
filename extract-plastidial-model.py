#!/usr/bin/env python3
import sys
import os
import argparse
from pathlib import Path

# ==============================================================================
# ModelSEED Database Dependency Check
# ==============================================================================
try:
    import BiochemPy 
except ImportError:
    print("\n" + "="*80)
    print(" CRITICAL ERROR: ModelSEEDDatabase Python Library Not Found!")
    print("="*80)
    print("This pipeline requires the 'ModelSEEDDatabase' repository to run.")
    print("\nINSTRUCTIONS:")
    print("1. Clone the repository: git clone https://github.com/ModelSEED/ModelSEEDDatabase.git")
    print("2. Export it: export PYTHONPATH=$PYTHONPATH:/absolute/path/to/ModelSEEDDatabase/Libs/Python/")
    print("="*80 + "\n")
    sys.exit(1)

project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/plastidial_model_extraction'
if module not in sys.path:
    sys.path.append(module)
from src.plastidial_model_extraction.plastidial_model_generator import ModelBuilder, ModelGenerator

def main():
    parser = argparse.ArgumentParser(description="Extracts a plastidial model from a full plant reconstruction.")
    parser.add_argument("-m", "--model", required=True, help="Path to the input full reconstruction file.")
    parser.add_argument("-o", "--output", help="Optional. Explicit output path. (Overrides the default 'plastidial-...' naming).")
    parser.add_argument("--media", default='data/metabolic_models/plastidial_biomass_media/PlantPlastidialAutotrophicMedia.json', help="Path to the media JSON.")
    parser.add_argument("--biomass", default='data/metabolic_models/plastidial_biomass_media/plastid_biomass.csv', help="Path to the biomass CSV.")
    
    args = parser.parse_args()

    print(f"Loading full model from: {args.model}")

    # 1. Initialize the Builder
    try:
        mBuilder = ModelBuilder(json_file=args.model)
    except TypeError:
        mBuilder = ModelBuilder(args.model)
        
    # 2. Intercept and Override the Output Name
    if args.output:
        out_path = args.output
        if out_path.endswith('.json'):
            out_path = out_path[:-5] 
        mBuilder.output_model_file = out_path
        print(f" -> Output destination overridden to: {args.output}")
    else:
        print(f" -> Using default output destination: {mBuilder.output_model_file}.json")

    # 3. Generate the core plastidial network
    mGen = ModelGenerator(mBuilder)
    mGen.generate_model()
    
    # 4. Patch in the missing exchange boundaries and biomass objective
    print(f" -> Applying media constraints from: {args.media}")
    mBuilder.load_media_file(args.media)
    mBuilder.add_media_exchange_reactions("_c0")
    
    print(f" -> Applying biomass equation from: {args.biomass}")
    mBuilder.add_biomass_bio1(args.biomass)

    # 5. Save the fully patched model natively
    print(" -> Saving the completed plastidial network natively...")
    mGen.clean_write_model(clean_up=False) 
    
    print("\nExtraction Complete!")

if __name__ == "__main__":
    main()