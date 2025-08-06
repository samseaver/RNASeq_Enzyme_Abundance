import sys
import os
from pathlib import Path
project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/reaction_scores'
if module not in sys.path:
    sys.path.append(module)
module = str(project_root)+'/src/plastidial_model_extraction'
if module not in sys.path:
    sys.path.append(module)

# plastidial model 
from src.plastidial_model_extraction.plastidial_model_generator import ModelBuilder, ModelGenerator

# Reaction scores 
import src.reaction_scores.computeScoresAndPredictions as csp

# Parameters file
from src.util.parameters import Parameters_QPSI

generate_models = False
compute_scores = True


if generate_models:
    spcs = ['Athaliana', 'Poplar', 'Sorghum']
    for spc in spcs:
        mBuilder = ModelBuilder(spc)
        mGen = ModelGenerator(mBuilder)
        mGen.run_model_generator(clean_up=True)

# Compute reaction scores
if compute_scores:
    project_param = Parameters_QPSI # 
    csp.generate_reactionScores(project_param)





