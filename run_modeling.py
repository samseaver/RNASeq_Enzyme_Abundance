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

# Parameters class
# The parameters dictate where to find the files
# And how to link them to different species/conditions
from src.util.parameters import Parameters_QPSI
from src.util.parameters import Parameters_ColdResponse
project_param = Parameters_ColdResponse()

# This flag is for the sake of our publication
# As described in our work, we take a full reconstruction
# And we reduce it to a working reconstruction of a plastid
compute_scores = True

if project_param.generate_plastidial_models:
    spcs = ['Athaliana', 'Poplar', 'Sorghum']
    for spc in spcs:
        mBuilder = ModelBuilder(spc)
        mGen = ModelGenerator(mBuilder)
        mGen.run_model_generator(clean_up=True)

# Compute reaction scores
if compute_scores:
    csp.generate_reactionScores(project_param,project_species=['TSU','C24'])
