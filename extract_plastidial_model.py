import sys
import os
from pathlib import Path

# ==============================================================================
# ModelSEED Database Dependency Check
# ==============================================================================
try:
    # Attempt to import the specific library you need from the database
    import BiochemPy 
except ImportError:
    print("\n" + "="*80)
    print(" CRITICAL ERROR: ModelSEEDDatabase Python Library Not Found!")
    print("="*80)
    print("This pipeline requires the 'ModelSEEDDatabase' repository to run.")
    print("Because this database is a data repository and cannot be installed via pip/conda,")
    print("you must clone it manually and point this script to it.")
    print("\nINSTRUCTIONS:")
    print("1. Clone the repository anywhere on your system:")
    print("   git clone https://github.com/ModelSEED/ModelSEEDDatabase.git")
    print("\n2. Point Python to the library folder. You can do this in two ways:")
    print("\n   Option A (In Code): Add these lines to the top of this script:")
    print("       import sys")
    print("       sys.path.append('/absolute/path/to/ModelSEEDDatabase/Libs/Python/')")
    print("\n   Option B (Terminal - Recommended): Export it to your PYTHONPATH:")
    print("       export PYTHONPATH=$PYTHONPATH:/absolute/path/to/ModelSEEDDatabase/Libs/Python/")
    print("="*80 + "\n")
    sys.exit(1)

project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/plastidial_model_extraction'
if module not in sys.path:
    sys.path.append(module)
from src.plastidial_model_extraction.plastidial_model_generator import ModelBuilder, ModelGenerator

# Load full model and generate plastidial model by extracting plastidial reactions
model_files = ['projects/qpsi/inputs/Sbicolor-v5.1-reconstruction_fixed.json', # Sorghum
               'projects/qpsi/inputs/Ptrichocarpa-v4.1-reconstruction_fixed.json'] # Poplar

media_file = 'data/metabolic_models/plastidial_biomass_media/PlantPlastidialAutotrophicMedia.json'
biomass_file = 'data/metabolic_models/plastidial_biomass_media/plastid_biomass.csv'

for model_file in model_files:
    mBuilder = ModelBuilder(json_file=model_file)
    mGen = ModelGenerator(mBuilder)
    mGen.run_model_generator(toJSON=False)

    # Add media exchange reactions from ModelSEED/KBase Media JSON
    mBuilder.load_media_file(media_file)
    mBuilder.add_media_exchange_reactions("_c0")

    # Add plastidial biomass from file
    mBuilder.add_biomass_bio1(biomass_file)

    # Write plastidial model to file
    mGen.clean_write_model(False)