import sys
import os
from pathlib import Path
project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/plastidial_model_extraction'
if module not in sys.path:
    sys.path.append(module)

# plastidial model 
from src.plastidial_model_extraction.plastidial_model_generator import ModelBuilder, ModelGenerator
spcs = ['Athaliana', 'Poplar', 'Sorghum']
for spc in spcs:
    mBuilder = ModelBuilder(spc)
    mGen = ModelGenerator(mBuilder)
    mGen.run_model_generator(clean_up=True)
