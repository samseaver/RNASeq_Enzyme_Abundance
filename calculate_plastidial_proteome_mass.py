import sys
import os
from pathlib import Path

project_root = Path(__file__).resolve().parent
module = str(project_root)+'/src/util'
if module not in sys.path:
    sys.path.append(module)

from src.util.proteinWeightGenerator import ProteinWeightGenerator

# compute protein weights
# pwg = ProteinWeightGenerator(spc.name, csp.project, spc.RNASeq_file_path,
#                             set(spc.metModel.modelfeatures_dict.keys()),
#                             csp.group_columns, cap_percent=percent)

totalProteinMass_df = pwg.plastid_weight_sums  # totalPlastidProteinMass()
if verbose: print(totalProteinMass_df.reset_index().head(5))
