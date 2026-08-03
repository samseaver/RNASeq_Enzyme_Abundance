import sys
import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
import csv
import json
from urllib.request import urlopen
import pandas as pa
import numpy as np
from scipy.stats import zscore
import math

import reactionScoresHelper as rsh
import fluxes_to_reactions as f2r

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from src.util.modelComponents import *

from src.util.proteinWeightGenerator import ProteinWeightGenerator
from parameters import *

import plotly.express as px
import matplotlib.pyplot as plt
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Model species information
class Species:
    def __init__(self, name, synonyms, modelJSON_file_path='', RNASeq_file_path='', model=None):
        self.name                = name
        self.synonyms            = synonyms
        self.metModel            = model
        self.pred_file_path      = ''
        self.modelJSON_file_path = modelJSON_file_path
        self.RNASeq_file_path    = RNASeq_file_path

    @property
    def modelJSON_file_path(self):
        return self._modelJSON_file_path

    @modelJSON_file_path.setter
    def modelJSON_file_path(self, modelJSON_file_path):
        self._modelJSON_file_path = modelJSON_file_path

        with open(modelJSON_file_path, 'r') as f:
            data=f.read()
            json_model = json.loads(data)
            self.metModel = Model().fromJSON(json_model)

            with open(self.name+'.json','w') as f:
                json.dump(self.metModel.as_dict(),f,indent=2)

        if self.metModel == None:
            raise ValueError(f"  Couldn't load metabolic model from {modelJSON_file_path}.")

## Pipeline parameters
class ComputeScoresPredictions:
    def __init__(self, param, project_species=list()):
        if param.error:
            sys.exit("Please fix errors above.")

        self.project = param.project
        self.project_species = project_species if project_species else param.project_species

        # TMM outlier capping removed: the percentile cap operated within the
        # metabolic-enzyme set and clipped the most abundant enzymes (Rubisco,
        # PPDK, carbonic anhydrase); transcript values are now used uncapped.

        self.PS_json_data = list()
        self.species_list =  list()
        self.cons_dict_list = list()

        self.value_column  = param.value_column
        self.modelJSON_files_folder = param.json_files_folder
        self.RNASeq_folder = param.RNASeq_folder
        # Write the results to this folder
        self.results_folder = param.results_folder
        self.reaction_file_suffix = param.reaction_file_suffix

        self.model_compartments = set()
        self.group_columns = param.group_columns
        self.rnaSeq_id_col = param.rnaSeq_id_col
        self.control_id = param.control_id
        self.trmt_colmn = param.trmt_colmn
        self.treatments = list()

# Read RNASeq data (transcript values used uncapped)
# Compute protein weights and total Plastid Protein Mass
# Compute relative molecular abundance (rma) to be used for relative-RES
def readRNASeq(spc, csp, verbose=True):
    print(f" Computing Protein weights for {spc.name} ...")
    # Read RANSeq data from file
    relab_df = pa.read_csv(spc.RNASeq_file_path)

    # # keep only genes involved in plastidial reactions
    relab_df = relab_df[(relab_df[csp.rnaSeq_id_col].isin(spc.metModel.modelfeatures_dict))]

    print("WE ARE READING RNASEQ FILE: ",spc.RNASeq_file_path)

    # TMM outlier capping removed --- transcript values are used uncapped.

    # compute protein weights
    pwg = ProteinWeightGenerator(spc.name, csp.project, spc.RNASeq_file_path,
                                set(spc.metModel.modelfeatures_dict.keys()),
                                csp.group_columns, cap_percent=-1)

    ## compute plastid protein mass and protein molecular weight
    pwg.computeWeightsAndMass(to_file=False)
    weights_df = pa.DataFrame.from_dict(pwg.protein_weight_dict, orient="index", columns=['weight'])
    weights_df = weights_df.rename_axis(csp.rnaSeq_id_col).reset_index()

    # Add gene weights
    relab_df = pa.merge(relab_df, weights_df, on=csp.rnaSeq_id_col, how='left')
    # get total plastid protein mass datafame
    totalProteinMass_df = pwg.plastid_weight_sums  # totalPlastidProteinMass()
    if verbose: print(totalProteinMass_df.reset_index().head(5))

    # Compute relative molecular abundance
    temp = pa.DataFrame()
    for name, group in relab_df.groupby(csp.group_columns):
        totalMass = totalProteinMass_df.loc[name[0], name[1], name[2]]['TotalPlastidMass']
        group['rma'] = (group['weight']*group[csp.value_column]) / totalMass
        group['rma_std'] = (group['weight']*group[csp.value_column+'_std']) / totalMass
        temp = pa.concat([temp, group])
    relab_df = temp

    temp = relab_df[["Gene_ID", "time_stamp", "value", "weight", "rma", "tissue", "treatment"]]
    temp = temp[(temp["tissue"] == 'Leaf') & (relab_df["treatment"] == 'Control')]
    temp["fold"] = temp["rma"] / temp["value"]
    ind = ['Gene_ID']
    col = ['time_stamp']
    val = ["value", "weight", "rma", "fold"]
    temp = temp.pivot(index=ind, columns=col, values=val)
    temp.columns = temp.columns.get_level_values(0) + '_' +  temp.columns.get_level_values(1)


    if verbose: print(relab_df.head(5))
    csp.value_column = 'rma'

    # #Compute molecular abundance: ralab / weight
    # relab_df['molab'] = relab_df['relab'] / relab_df['weight']
    # for name, group in relab_df.groupby(csp.group_columns):
    #     totalMass = totalProteinMass_df.loc[name[0], name[1], name[2]]['TotalPlastidMass']
    #     relab_df['relab'] = relab_df['relab'] * totalMass

    # Save dataframe to file
    relab_path = os.path.join(csp.results_folder, f"{spc.name}_relab.csv")
    relab_df.to_csv(relab_path, index=False)

    if verbose: print(relab_df.head(5))
    return relab_df, pwg

# Read Trimmed Mean of the M-values (TMM) values from file.
# Read TMM values (used uncapped)
def readTMMdata(spc, csp, verbose=True):
    # Read RANSeq data from file
    tmm_file = spc.RNASeq_file_path

    # TSV format is canonical; pandas auto-decompresses .gz/.xz/.bz2 by extension
    tmm_df = pa.read_csv(tmm_file, sep='\t')

    # Check for ID of first column for gene/transcript/protein ids
    # Sometimes it is empty or has a different name, best to keep it consistent
    if(tmm_df.columns[0] != csp.rnaSeq_id_col):
        tmm_df.columns.values[0] = csp.rnaSeq_id_col
        
    # keep only genes involved in plastidial reactions
    if verbose:
        print(f"\tThere are {len(spc.metModel.modelfeatures_dict)} genes in the model")
        print(f"\t          {len(tmm_df[csp.rnaSeq_id_col].unique())} genes in the RNAseq")
    tmm_df = tmm_df[(tmm_df[csp.rnaSeq_id_col].isin(spc.metModel.modelfeatures_dict))]

    # TMM outlier capping removed --- values are used uncapped.

    return tmm_df

def generate_reactionScores(parameters, project_species:list=[], verbose=False):
    csp = ComputeScoresPredictions(parameters, project_species)

    method = 'sum'
    print("Using default Sum-Min-Sum approach to compute reaction scores")
    
    # set the file paths to RNAseq and models
    json_dict = parameters.model_paths
    if verbose:
        for species, json_file in json_dict.items():
            print("        "+species+" <-- "+json_file)
    rnaseq_dict = parameters.rnaseq_paths
    if verbose:
        for species, rnaseq_file in rnaseq_dict.items():
            print("        "+species+" <-- "+rnaseq_file)
            
    # set model species
    csp.species_list = [Species(spc, parameters_all_spc[spc]['synonyms'], json_dict[spc], rnaseq_dict[spc]) for spc in csp.project_species]

    for species in csp.species_list:

        if verbose:
            print("\n"+"-"*10+"    Computing reaction scores for "+species.name+"-"*10)
        ## Start reaction score computation using TMM -- used for K_app computation -------------
        if verbose:
            print("    Computing reactions scores using TMM data")
        tmm_df = readTMMdata(species, csp,verbose=verbose)

        if 'value_log' not in tmm_df.columns:
            tmm_df['value_log'] = np.log(tmm_df[csp.value_column])

        ## Compute Reaction Scores ============================================================================
        reaction_scores = rsh.compute_model_score(tmm_df, parameters, species, csp, method=method, verbose=verbose)
        ##=====================================================================================================

        if 'reaction_id' not in reaction_scores.columns:
            print(f"Error: Zero reaction scores were generated.")
            print(f"This is an indication that none of the genes in the model were linked to transcripts in the RNASeq dataset")
            sys.exit(1)
    
        # --- Statistics Reporting ---
        total_model_rxns = len(species.metModel.modelreactions_dict)
        
        # Count reactions that actually have gene mappings (GPRs)
        gene_associated_rxns = 0
        for rxn in species.metModel.modelreactions_dict.values():
            if rxn.genes: # Checks if the gene list is not empty
                gene_associated_rxns += 1
                
        # Count reactions that successfully received a score
        scored_unique_rxns = reaction_scores['reaction_id'].nunique()
        
        print(f"\n" + "="*50)
        print(f"  Reaction Score Stats: {species.name}")
        print(f"  1. Total Model Reactions:       {total_model_rxns}")
        print(f"  2. Gene-Associated Reactions:   {gene_associated_rxns}")
        print(f"  3. Reactions with Transcript:   {scored_unique_rxns}")
        
        if total_model_rxns > 0:
            pct_gene = (gene_associated_rxns / total_model_rxns) * 100
            print(f"     -> Gene Coverage of Model:   {pct_gene:.2f}%")
            
        if gene_associated_rxns > 0:
            pct_score = (scored_unique_rxns / gene_associated_rxns) * 100
            print(f"     -> Transcript Coverage:      {pct_score:.2f}% (of gene-associated)")
            
        print("="*50 + "\n")
        # ----------------------------

        # RENAME columns for clarity
        reaction_scores.rename(columns={
            csp.value_column: 'reaction_score',    # Was 'mean_value'
            csp.rnaSeq_id_col: 'limiting_subunit'  # Was 'Geneid'
        }, inplace=True)

        # Reorder columns and implicitly drop anything not in this list (like 'features')
        reaction_scores = reaction_scores[['condition', 'reaction_id', 'reaction_score', 'limiting_subunit']]

        #  Save results to file
        csp.objSaveTo = os.path.join(csp.results_folder,f"{species.name}{csp.reaction_file_suffix}")
        reaction_scores.to_csv(csp.objSaveTo, index=False, sep='\t')
        print("Reaction scores saved to",csp.objSaveTo)

if __name__ == '__main__':
    #Parameters_test_a, Parameters_QPSI, Parameters_secMeta_TSU
    project = Parameters()
    generate_reactionScores(project)

