from pathlib import Path
import os
import sys

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
runtime_root = os.getcwd()
from src.util.bcolors import bcolors

# Even though the idea is to sort out data by species
# If you have data from different genotypes/ecotypes
# Of the same species, you will need to separate those out too
# If they would of course use the same full species reconstruction
# You can put the species in the synonyms
parameters_all_spc={"TSU":{'name':'TSU','synonyms':['TSU','Athaliana']},
                    "C24":{'name':'C24','synonyms':['C24','Athaliana']},
                    "Athaliana":{'name': "Arabidopsis",
                                 'synonyms': ["Atha", "athaliana", "Arabidopsis", 'sandbox','Athaliana']},
                    "Poplar":{'name': "Poplar",
                              'synonyms': ["Poplar", "Pptrich", "Ptrichocarpa", "ptrich", "ptr"]},
                    "Sorghum":{'name': "Sorghum",
                               'synonyms': ["Sorghum", "Sbicolor","Sbi","sbi"]}}

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

        self.objSaveTo       = ''
        self.relabSaveTo     = ''

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

# Sets project specif information for automatic processing of transcriptome
class Parameters_ColdResponse:
    def __init__(self):
        self.project         = 'ColdResponse'
        # This parameter allows you to switch between different value types
        # But they should be saved in a folder of the same name in the RNASeq_folder
        # i.e. .../rnaseq-data/tmm/
        self.msr             = 'tmm' # 'counts', 'tmm', 'DESeq2' or 'tpm'
        self.value_column    = 'mean_value'
        self.std_column      = 'std_value'
        self.timepoints      = []
        self.project_species = ["Atha"]
        self.group_columns   = ['treatment', 'time_stamp']

        # This is the name of the column that contains gene/transcript/protein ids
        self.rnaSeq_id_col   = 'gene_id'
        self.error           = False
        self.control_id      = 'CTL'
        self.trmt_colmn      = 'treatment'

        self.objSaveTo       = ''
        self.relabSaveTo     = ''

        self.generate_plastidial_models = False
        
        #Load the models from here
        self.json_files_folder = runtime_root+"/projects/cold-response/inputs/"
        if not os.path.exists(self.json_files_folder):
            print("Please move all species models to:")
            print(self.json_files_folder)
            self.error = True

        # RNAseq data will be here
        self.RNASeq_folder = runtime_root+"/projects/cold-response/rnaseq-data/"
        if not os.path.exists(self.RNASeq_folder):
            print("Please move RNASeq data to:")
            print(self.RNASeq_folder)
            self.error = True

        # Write the results to this folder
        self.results_folder = runtime_root+"/projects/cold-response/integration-results/"
        if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)

        # Access PlantSEED DB
        self.PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
        self.PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
        self.PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"

class Parameters_BRaVE:
    def __init__(self):
        self.project         = 'brave'
        # This parameter allows you to switch between different value types
        # But they should be saved in a folder of the same name in the RNASeq_folder
        # i.e. .../rnaseq-data/tmm/
        self.msr             = 'tmm' # 'counts', 'tmm', 'DESeq2' or 'tpm'
        self.value_column    = 'mean_value'
        self.std_column      = 'std_value'
        self.timepoints      = []
        self.project_species = ["Sbicolor"]
        self.group_columns   = ['treatment']

        # This is the name of the column that contains gene/transcript/protein ids
        self.rnaSeq_id_col   = 'Geneid'
        self.error           = False
        self.control_id      = 'CTL'
        self.trmt_colmn      = 'treatment'

        self.objSaveTo       = ''
        self.relabSaveTo     = ''

        self.generate_plastidial_models = False
        
        #Load the models from here
        self.json_files_folder = runtime_root+f"/projects/{self.project}/inputs/"
        if not os.path.exists(self.json_files_folder):
            print("Please move all species models to:")
            print(self.json_files_folder)
            self.error = True

        # RNAseq data will be here
        self.RNASeq_folder = runtime_root+f"/projects/{self.project}/rnaseq-data/"
        if not os.path.exists(self.RNASeq_folder):
            print("Please move RNASeq data to:")
            print(self.RNASeq_folder)
            self.error = True

        # Write the results to this folder
        self.results_folder = runtime_root+f"/projects/{self.project}/integration-results/"
        if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)

        # Access PlantSEED DB
        self.PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
        self.PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
        self.PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"
