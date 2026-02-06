import sys
import os
from pathlib import Path
project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/reaction_scores'
if module not in sys.path:
    sys.path.append(module)

# Parameters class
# The parameters dictate where to find the files
# And how to link them to different species/conditions
from src.util.parameters import Parameters
project_param = Parameters()

# Reaction scores 
import src.reaction_scores.computeScoresAndPredictions as csp
for species in project_param.project_species:
    csp.generate_reactionScores(project_param,project_species=[species])
