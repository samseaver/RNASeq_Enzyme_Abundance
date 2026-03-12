from pathlib import Path
import os
import sys

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
runtime_root = os.getcwd()

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
class Parameters:
    def __init__(self):
        self.project         = 'qpsi-plastidial'
        self.msr             = 'tmm'
        self.value_column    = 'value'
        self.project_species = ["Sorghum","Poplar"]
        self.group_columns   = ['condition']
        self.rnaSeq_id_col   = 'Gene_ID'
        self.error           = False
        self.control_id      = 'Control'
        self.trmt_colmn      = 'condition'

        # if method = 'sum' then compute reaction scores using the sum-min-sum method
        # if method = 'relab' then compute reaction scores in the same manner
        #   but normalized to the plastidial proteome
        self.reaction_score_method = 'sum'

        #Load the models from here
        self.json_files_folder = project_root+"projects/"+self.project+"/inputs/"
        if not os.path.exists(self.json_files_folder):
            print("Please move all species models to:")
            print(self.json_files_folder)
            self.error = True

        # RNAseq data will be here
        self.RNASeq_folder = project_root+"projects/"+self.project+"/rnaseq-data/"
        if not os.path.exists(self.RNASeq_folder):
            print("Please move RNASeq data to:")
            print(self.RNASeq_folder)
            self.error = True

        # Write the results to this folder
        self.results_folder = project_root+"projects/"+self.project+"/integration_results/"
        if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)

        self.ignore_organellar_roles = 'data/organellar-encoded_subunits_to_ignore.txt'
        
        # Access PlantSEED DB
        self.PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
        self.PS_tag  = "dev"
        self.PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"
