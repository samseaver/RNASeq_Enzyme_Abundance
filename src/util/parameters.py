from pathlib import Path
import os
import sys

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)

from src.util.bcolors import bcolors


parameters_all_spc={"Atha":
                        {'name': "Atha",
                                'synonyms': ["Atha", "athaliana", "Arabidopsis", 'sandbox']},
                    "Poplar":
                        {'name': "Poplar",
                                'synonyms': ["Poplar", "Pptrich", "Ptrichocarpa", "ptrich", "ptr"]},
                    "Sorghum":
                        {'name': "Sorghum",
                                'synonyms': ["Sorghum", "Sbicolor_"]},
                    
                    # Add species of interest here 
                    }

# Sets project specif information for automatic processing of transcriptome
class Parameters_QPSI:
    def __init__(self):
        self.project         = 'QPSI'
        self.msr             = 'tmm' # 'counts', 'tmm', 'DESeq2' or 'tpm'
        self.value_column    = 'molab' if (self.msr == 'counts') else 'value' # 'relab'
        self.metals          = ['Fe', 'Zn']
        self.project_species = ["Atha", "Sorghum", "Poplar"]
        self.group_columns   = ['tissue', 'treatment', 'time_stamp']
        self.rnaSeq_id_col   = 'Gene_ID'
        self.error           = False
        self.control_id      = 'Control'
        self.trmt_colmn      = 'treatment'

        project_root = str(Path(__file__).resolve()).split('src')[0]

        # # Get predictions from here
        # self.rgb_folders = project_root+"data/prediction_folders"
        # if not os.path.exists(self.rgb_folders):
        #     print("Please move prediction folders for all species and metals to:")
        #     print(bcolors.FAIL+self.rgb_folders+bcolors.ENDC)
        #     print(bcolors.WARNING+" Note: Folders are expected to have the following structure:\n\
        #     name_specifying_<species>_<metal>/<species>/[all domains]"+bcolors.ENDC)
        #     self.error = True

        #Load the models from here
        self.json_files_folder = project_root+"data/metabolic_models/plastidial_models"
        # self.json_files_folder = project_root+"data/metabolic_models/fullmodels_media"
        if not os.path.exists(self.json_files_folder):
            print("Please move all species models to:")
            print(bcolors.FAIL+self.json_files_folder+bcolors.ENDC)
            self.error = True

        # RNAseq data will be here
        self.RNASeq_folder = project_root+"data/RNASeq_data/"
        if not os.path.exists(self.RNASeq_folder):
            print("Please move RNASeq data to:")
            print(bcolors.FAIL+self.RNASeq_folder+bcolors.ENDC)
            self.error = True

        # Write the results to this folder
        self.results_folder = project_root+"integration_results/"
        if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)

        # Access PlantSEED DB
        self.PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
        self.PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
        self.PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"

