import sys
import warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
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
# from src.metalBinding.metalRdgbHelper import RDGBHelper
from src.util.proteinWeightGenerator import ProteinWeightGenerator
from src.util.bcolors import bcolors as bc
from src.util.parameters import *

import plotly.express as px
import matplotlib.pyplot as plt
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#Global variables
PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"

def generateRolesDict():
    print(bc.SUBRESULT+"PS Roles reading ...")
    # PS_json_data = json.load(urlopen(PS_url+PS_tag+PS_json))

    roles_file = os.path.join(project_root, "data", "metabolic_models", "PlantSEED_Roles.json")
    with open(roles_file, 'r') as f:
        PS_json_data = json.load(f)

    roles_dict = dict()
    for item in PS_json_data:
        # Create a new role entry for the consolidated model
        role_dict = dict()
        role_dict["role"]        = item["role"]
        role_dict["subsystems"]  = item["subsystems"]
        role_dict["reactions"]   = list()
        roles_dict[item["role"]] = role_dict

    i = 1
    for fr in ["1", "2", "3", "C1", "C2"]:
        role_dict = dict()
        role_dict["role"]        = f"Ferredoxin {fr} (EC 0.0.0.{i})"
        role_dict["subsystems"]  = ["Ferredoxins"]
        role_dict["reactions"]   = [f"rxn00fd{fr}_d0", "rxn00fdA_d0"]
        roles_dict[role_dict["role"]] = role_dict
        i += 1
    return roles_dict
roles_dict = generateRolesDict()


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
            self.metModel = Model().fromJSON(json_model, roles_dict)

        if self.metModel == None:
            raise ValueError(bc.FAIL+f"  Couldn't load metabolic model from {modelJSON_file_path}."+bc.ENDC)

    def updateMetalBinding(self, id_pfam_dict, metal):
        intersection = set(id_pfam_dict.keys()).\
                        intersection(self.metModel.modelfeatures_dict.keys())

        for gene_id in intersection:
            self.metModel.modelfeatures_dict[gene_id].binds = metal
            self.metModel.modelfeatures_dict[gene_id].pfam.add(id_pfam_dict[gene_id])


## Pipeline parameters
class ComputeScoresPredictions:
    def __init__(self, parameters_class, project='QPSI', project_species=list()):
        param = parameters_class()

        if param.error:
            sys.exit("Please fix errors above.")

        self.project = param.project
        self.project_species = project_species if project_species else param.project_species
        self.outlier_cap_percentile = 95

        # if not os.path.exists(param.rgb_folders):
        #     raise FileNotFoundError("Prediction folders not found.")

        # self.folders = os.listdir(param.rgb_folders)
        # self.folders = [os.path.join(param.rgb_folders, f) for f in self.folders]

        self.PS_json_data = list()
        self.species_list =  list()
        self.cons_dict_list = list()

        self.msr = param.msr
        self.metals = param.metals
        self.value_column  = param.value_column
        self.modelJSON_files_folder = param.json_files_folder
        self.RNASeq_folder = param.RNASeq_folder
        # Write the results to this folder
        self.results_folder = param.results_folder

        self.model_compartments = set()
        self.group_columns = param.group_columns
        self.rnaSeq_id_col = param.rnaSeq_id_col
        self.control_id = param.control_id
        self.trmt_colmn = param.trmt_colmn
        self.treatments = list()

        self.objSaveTo = param.objSaveTo
        self.relabSaveTo = param.relabSaveTo


# Find the path to JSON model files for each species
def jsonModelPath2Var(csp, verbose=True):
    json_dict = dict()
    if verbose: print(bc.PROMPT+"-"*10+" Updating JSON file paths "+"-"*10)
    fileNames = [os.path.join(csp.modelJSON_files_folder, x)
                    for x in os.listdir(csp.modelJSON_files_folder)]

    for fileName in fileNames:
        for spc in csp.project_species:
            synonyms = parameters_all_spc[spc]['synonyms']
            if any(y.lower() in fileName.lower() for y in synonyms) and fileName.endswith('.json'):
            #  \
            # and (('070224'.lower() in fileName.lower()) or ('250512'.lower() in fileName.lower())):
            # Athaliana_Thylakoid_Reconstruction_ComplexFix_070224
                json_dict[spc] = fileName
                break

    return json_dict


# Find the path to RNASeq files for each species
def rnaSeqFile2Var(csp, verbose=True):
    rnaseq_dict = dict()

    if verbose: print(bc.PROMPT+"-"*10+" Updating RNASeq file paths "+"-"*10)
    fileNames = os.listdir(csp.RNASeq_folder+csp.msr+"/")
    for fileName in fileNames:
        if (csp.msr in fileName.lower()) and ('mean' in fileName.lower()):# and ('counts' in fileName.lower()):
            for spc in csp.project_species:
                if (csp.project != 'secMeta') and ('atha' in spc.lower()):
                    rnaseq_dict[spc] = ""
                    continue

                synonyms = parameters_all_spc[spc]['synonyms']
                if any(y.lower() in fileName.lower() for y in synonyms):
                    rnaseq_dict[spc] = csp.RNASeq_folder+csp.msr+"/"+fileName
                    break

    return rnaseq_dict

# Read RNASeq data, remove outliers by caping at the "outlier_cap_percentile"
# Compute protein weights and total Plastid Protein Mass
# Compute relative molecular abundance (rma) to be used for relative-RES
def readRNASeq(spc, csp, verbose=True):
    print(bc.PROG+f" Computing Protein weights for {spc.name} ..."+bc.ENDC)
    # Read RANSeq data from file
    relab_df = pa.read_csv(spc.RNASeq_file_path)

    # # keep only genes involved in plastidial reactions
    relab_df = relab_df[(relab_df[csp.rnaSeq_id_col].isin(spc.metModel.modelfeatures_dict))]

    # format time_stamp
    if csp.project == 'QPSI':
        # add '0' to single digit time points for ordering purposes
        relab_df['time_stamp'] = relab_df['time_stamp'].transform(lambda ts:ts if ts in ['14d', '21d'] else '0'+ts)


    # # Remove outliers: !!! NOT NEEDED WHEN USING COUNTS (ALREADY CAPPED) !!!
    percent = -1
    if csp.msr in ['tpm', 'tmm']:
        percent = relab_df[csp.value_column].\
            describe([csp.outlier_cap_percentile/100])\
            [str(csp.outlier_cap_percentile)+'%']

        relab_df.loc[(relab_df[csp.value_column] > percent), csp.value_column] = percent

    # compute protein weights
    pwg = ProteinWeightGenerator(spc.name, csp.project, spc.RNASeq_file_path,
                                set(spc.metModel.modelfeatures_dict.keys()),
                                csp.group_columns, cap_percent=percent)

    ## compute plastid protein mass and protein molecular weight
    pwg.computeWeightsAndMass(to_file=False)
    weights_df = pa.DataFrame.from_dict(pwg.protein_weight_dict, orient="index", columns=['weight'])
    weights_df = weights_df.rename_axis(csp.rnaSeq_id_col).reset_index()
    print(weights_df.head())
    # print(abc)

    # Add gene weights
    relab_df = pa.merge(relab_df, weights_df, on=csp.rnaSeq_id_col, how='left')
    # get total plastid protein mass datafame
    totalProteinMass_df = pwg.plastid_weight_sums  # totalPlastidProteinMass()
    if verbose: print(totalProteinMass_df.reset_index().head(5))
    print(totalProteinMass_df.reset_index().head(5))

    print(abc)
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

    # temp.to_csv("compare_weighted_value_control.csv")
    # print(abc)

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
# Cap the values at "outlier_cap_percentile" to remove outliers
def readTMMdata(spc, csp, verbose=True):
    # Read RANSeq data from file
    tmm_file = spc.RNASeq_file_path
    if csp.project == 'QPSI':
        tmm_file = tmm_file.replace(csp.msr, 'tmm')
    tmm_df = pa.read_csv(tmm_file)


    # keep only genes involved in plastidial reactions
    print(bc.PROG +f"\tThere are {len(spc.metModel.modelfeatures_dict)} genes in the model")
    print(bc.PROG +f"\t          {len(tmm_df[csp.rnaSeq_id_col].unique())} genes in the RNAseq")
    tmm_df = tmm_df[(tmm_df[csp.rnaSeq_id_col].isin(spc.metModel.modelfeatures_dict))]

    # Cap TMM values at the outlier_cap_percentile percentile
    percent = tmm_df[csp.value_column].\
        describe([csp.outlier_cap_percentile/100])\
        [str(csp.outlier_cap_percentile)+'%']
    tmm_df.loc[(tmm_df[csp.value_column] > percent), csp.value_column] = percent

    # Histogram and violin plot to double check
    bin_width= 10
    wt, ht = 800, 600
    fig = px.violin(tmm_df,
                    x='time_stamp',
                    y='value',
                    color="treatment",
                    title=spc.name+" after cap",
                    height=ht, width=wt
                )
    # fig.update_yaxes(range=[-1500, 80000])
    fig.show()

    nbins = math.ceil((tmm_df[csp.value_column].max() - tmm_df[csp.value_column].min()) / bin_width)
    fig = px.histogram(tmm_df, x=csp.value_column, nbins=nbins)
    fig.show()


    # Format time points 
    if csp.project == 'QPSI':
        # add '0' to single digit time points
        tmm_df['time_stamp'] = tmm_df['time_stamp'].transform(lambda ts:ts if ts in ['14d', '21d'] else '0'+ts)


    ## -----------------------------------------------------------------------------------------
    if 'atha' in spc.name.lower():
        # remove ATC and ATM
        if verbose: print('-----> before removing the ATC/ATM ', tmm_df.shape)
        clm = 'Gene_Id' if 'Gene_Id' in tmm_df else 'Gene_ID'
        filter = tmm_df[clm].str.contains("ATC.*|ATM.*")
        tmm_df = tmm_df[~filter]
        if verbose: 
            print('-----> after removing the ATC/ATM ', tmm_df.shape)
            print(tmm_df.head(5))

    return tmm_df

def generate_reactionScores(parameters_class, project:str='QPSI', project_species:list=[], verbose=True):
    csp = ComputeScoresPredictions(parameters_class, project, project_species)

    if not os.path.exists(csp.results_folder):
        os.makedirs(csp.results_folder)

    # set the file paths to RNAseq and models
    json_dict = jsonModelPath2Var(csp)
    if verbose:
        for spc, json_file in json_dict.items():
            print(bc.PROG+"        "+spc+" <-- "+json_file+bc.ENDC)
    rnaseq_dict = rnaSeqFile2Var(csp)
    if verbose:
        for spc, rnaseq_file in rnaseq_dict.items():
            print(bc.PROG+"        "+spc+" <-- "+rnaseq_file+bc.ENDC)
    # set model species
    csp.species_list = [Species(spc, parameters_all_spc[spc]['synonyms'], json_dict[spc], rnaseq_dict[spc]) for spc in csp.project_species]


    ## Reactions cores ============================================================================
    for spc in csp.species_list:

        print(bc.PROMPT+"\n"+"-"*10+"    Computing reaction scores "+"-"*10+bc.ENDC)
        if ('atha' in spc.name.lower()) and (csp.project != 'secMeta'): continue
        ## Start reaction score computation using TMM -- used for K_app computation -------------
        print(bc.PROG +"    Computing reactions scores using TMM data"+bc.ENDC)
        tmm_df = readTMMdata(spc, csp)
        if 'value_log' not in tmm_df.columns:
            tmm_df['value_log'] = np.log(tmm_df[csp.value_column]) # np.log(tmm_df['value'])

        csp.treatments = list(tmm_df[csp.trmt_colmn].unique())
        csp.treatments.remove(csp.control_id)
        #  Compute reaction scores using the sum-min-sum method
        tmm_scores = rsh.compute_model_score(tmm_df, spc.metModel, csp.rnaSeq_id_col,\
                                 csp.value_column, csp.group_columns, spc.name, method='sum', verbose=True)

        tmm_scores = rsh.compute_rxn_variability(tmm_scores, csp.group_columns, csp.treatments, csp.control_id, csp.trmt_colmn, value_col=csp.value_column, percentile=90, verbose=False)

        #  Save results to file
        if not csp.objSaveTo: 
            csp.objSaveTo = f"{csp.results_folder}{spc.name}_objective_abundance_{csp.control_id}.tsv"
        tmm_scores.to_csv(csp.objSaveTo, index=False, sep='\t')
        # continue

        if csp.project not in ["QPSI", "hAlpha"]:
            sys.exit(f"Scores saved to: \n {csp.results_folder}{spc.name}_objective_abundance_{csp.control_id}.csv")
        ###### --- STOP HERE IF RELATIVE RS IS NOT RELEVANT 


        ## RNASeq data and protein weight processing --------------------------------------------
        print(bc.SUBRESULT +"        ** {} Model Reactions.".format(spc.name)+bc.ENDC)
        relab_df, pwg = readRNASeq(spc, csp)

        ## Start relative abundance based reaction score computation ----------------------------
        print(bc.PROG +"    Computing reactions scores"+bc.ENDC)
        # set treatment list from the new DF in case they are different from previously set values
        csp.treatments = list(relab_df[csp.trmt_colmn].unique())
        csp.treatments.remove(csp.control_id)

        # compute the relative abundance scores -------------------------------------------------
        rxn_scores = rsh.compute_model_score(relab_df, spc.metModel, csp.rnaSeq_id_col,\
                csp.value_column, csp.group_columns, spc.name, method='relab', verbose=True)

        # Add fluxes to each reaction -----------------------------------------------------------
        print(bc.PROMPT+"    Running COBRApy FVA"+bc.ENDC)
        fva_df = pa.DataFrame()
        merge_on = ['rxn_ID', 'tissue'] if csp.project == 'QPSI' else ['rxn_ID']

        #  Compute the min and max fluxes using Flux Variability Analysis (FVA)
        #  Use different models for root and leaf to generate FVA values for QPSI model
        for tissue in ['Leaf']:#, 'Root']:
            model_file = spc.modelJSON_file_path.replace('.json', '.xml')
            if tissue == 'Root':
                if csp.project != 'QPSI':
                    continue
                model_file = model_file.replace('.xml', '_sucrose.xml')

            fva_df = pa.concat([fva_df, f2r.run_FVA(model_file, tissue, spc.name)], ignore_index=True)

        if not fva_df.empty:
            if csp.project != 'QPSI':
                fva_df.drop(['tissue'], axis=1, inplace=True)
            rxn_scores =  pa.merge(rxn_scores, fva_df, how="left", on=merge_on)

        # Setting flexibility
        rxn_scores['flexibility'] = 'none'
        rxn_scores.loc[np.abs(rxn_scores['maximum']-rxn_scores['minimum']) <= 10**-6, 'flexibility'] = 'fixed'
        rxn_scores.loc[np.abs(rxn_scores['maximum']-rxn_scores['minimum']) > 10**-6, 'flexibility'] = 'flexible'

        # Add dist to identity line -------------------------------------------------------------
        csp.value_column =  "value"
        if "weightedSumWeight" in rxn_scores.columns:
            rxn_scores.rename(columns={'weightedSumWeight': csp.value_column}, inplace=True)

        print(rxn_scores.columns)
        rxn_scores = rsh.compute_rxn_variability(rxn_scores, csp.group_columns, csp.treatments, csp.control_id, csp.trmt_colmn, csp.value_column, percentile=90, verbose=False)

        print(bc.PROMPT+"Writing reaction scores and details to:")
        print(csp.results_folder+spc.name+"_rxn_scores_"+csp.msr+".csv"+bc.ENDC)

        if not csp.relabSaveTo:
            csp.relabSaveTo = csp.results_folder+spc.name+"_relab_rxn_scores_"+csp.msr+".csv"
        rxn_scores.to_csv(csp.relabSaveTo, index=False)

        #  Processing file genrates a new scores file that edits the format of the DF
        #   if the pipeline is run again, remove the existing edited RES file in order to
        #   compute a new one in the processing script.
        rxn_scores_edited = csp.results_folder+"{}_rxn_scores_{}_edited.csv".format(spc.name, csp.msr)
        if os.path.exists(rxn_scores_edited):
            os.remove(rxn_scores_edited)

        #  Write species list of features and their binding properties to be used for result
        #    processing
        with open(csp.results_folder+spc.name+"_all_features.csv", "w") as outfile:
            all_features = spc.metModel.modelfeatures_dict.keys()
            for key, ftr in spc.metModel.modelfeatures_dict.items():
                bind = None if not ftr.binds else "_".join(ftr.binds)
                outfile.write(f"{key}, {bind}\n")
        # continue
        # Compare RES computation methods ------------------------------------------------------
        plotMethodDiffs(csp, spc, tmm_df, tmm_scores, rxn_scores, pwg)

## Compare different reaction score computation methods
##  (sum_min_sum) vs. (max_min_max) vs. relative abundance
##  plot values as bar plots, scatter plots, and histogram
def plotMethodDiffs(csp, spc, tmm_df, tmm_scores, relative_scores, pwg, verbose=False):
    ## for figures ########################### ########################### ###############
    tmm_old_scores = rsh.compute_model_score(tmm_df, spc.metModel, csp.rnaSeq_id_col,\
        csp.value_column, csp.group_columns, spc.name, method='max', verbose=True)
    tmm_old_scores = rsh.compute_rxn_variability(tmm_old_scores, csp.group_columns, csp.treatments, csp.control_id, csp.trmt_colmn, value_col='value', percentile=90, verbose=False)

    tmm_old_scores = tmm_old_scores[['rxn_ID', 'value', 'rxn_score_I_dist',
                                'rxn_dist_quantile']+csp.group_columns]
    tmm_old_scores.rename(columns={'value': 'value_max',
                        'rxn_score_I_dist':'rxn_score_I_dist_max',
                        'rxn_dist_quantile': 'rxn_dist_quantile_max'},
                        inplace=True)

    tmm_old_scores =  pa.merge(tmm_old_scores, tmm_scores[['rxn_ID', 'value',
                                'rxn_score_I_dist', 'rxn_dist_quantile']+csp.group_columns],
                                how="inner", on=['rxn_ID']+csp.group_columns)
    tmm_old_scores.rename(columns={'value': 'value_sum',
                        'rxn_score_I_dist':'rxn_score_I_dist_sum',
                        'rxn_dist_quantile': 'rxn_dist_quantile_sum'},
                        inplace=True)

    tmm_old_scores =  pa.merge(tmm_old_scores, relative_scores[['rxn_ID', 'value','rxn_score_I_dist', 'rxn_dist_quantile']+csp.group_columns],
                                how="inner", on=['rxn_ID']+csp.group_columns)
    tmm_old_scores.rename(columns={'value': 'value_relative',
                        'rxn_score_I_dist':'rxn_score_I_dist_relative',
                        'rxn_dist_quantile': 'rxn_dist_quantile_relative'},
                        inplace=True)
    if verbose: print(tmm_old_scores.head(5))


    # scores_df['norm_'+value_col] =  scores_df.groupby(gcs)[value_col]\
    #     .apply(lambda value : (value - value.min()) / (value.max() - value.min()))
    tmm_old_scores = tmm_old_scores[tmm_old_scores['tissue'] == 'Leaf']
    temp = pa.DataFrame()
    for name, group in tmm_old_scores.groupby(['treatment', 'time_stamp']):
        mx = max(group['value_sum'].max(), group['value_max'].max())
        mn = min(group['value_sum'].min(), group['value_max'].min())
        group['value_sum'] = (group['value_sum'] - mn) / (mx - mn)
        group['value_max'] = (group['value_max'] - mn) / (mx - mn)

        temp =  pa.concat([group, temp], ignore_index=True)

    tmm_old_scores = temp
    # tmm_old_scores[['value_sum', 'value_max']] =  tmm_old_scores.groupby(csp.group_columns)[['value_sum', 'value_max']]\
    #     .apply(lambda value : (value - value.min()) / (value.max() - value.min()))

    # Create traces
    tmm_old_scores = tmm_old_scores[~(tmm_old_scores['time_stamp'] == '00h')]

    tmm_old_scores = pa.melt(tmm_old_scores, id_vars=['rxn_ID', 'treatment', 'time_stamp'], value_vars=['rxn_dist_quantile_sum', 'rxn_dist_quantile_max', 'rxn_dist_quantile_relative'])

    # $african-violet: 'rgba(171, 146, 191, 1)';
    # $hookers-green: 'rgba(78, 110, 93, 1)';
    # $apricot: 'rgba(255, 211, 186, 1)';
    # $blue-munsell: rgba(29, 138, 153, 1);
    # $jasmine: rgba(251, 216, 127, 1);
    # $light-orange: rgba(252, 208, 161, 1);
    color_map = {'rxn_dist_quantile_max': 'rgba(171, 146, 191, 1)', 'rxn_dist_quantile_sum': 'rgba(29, 138, 153, 1)', 'rxn_dist_quantile_relative': 'rgba(252, 208, 161, 1)'}
    method_map = {'rxn_dist_quantile_max': 'Reaction score w/ max', 
                  'rxn_dist_quantile_sum': 'Reaction score w/ sum', 
                  'rxn_dist_quantile_relative': 'Plastid-relative reaction score w/ sum'}

    method_map = {'rxn_dist_quantile_max': r"$r_{s} \; max$", 
                  'rxn_dist_quantile_sum': r"$r_{s} \; sum$", 
                  'rxn_dist_quantile_relative': r"$r_{s}^{r} \; sum$"}

    counts_df = tmm_old_scores[tmm_old_scores['value']>=80].groupby(['treatment', 'time_stamp', 'variable']).size().reset_index(name='counts')
    print("-------- ", spc.name)
    # print(counts_df[counts_df['treatment'].isin(['FeLim', 'Control'])])
    print(counts_df[counts_df['treatment'].isin(['FeLim', 'Control'])].groupby(['variable'])['counts'].sum())
    # print(abc)

    if verbose: print(counts_df.head(4))
    treatments = list(counts_df['treatment'].unique())
    treatments = ['FeLim']#, 'FeEX', 'ZnLim', 'ZnEx']
    species = [spc.name]
    specs = []
    for trm in treatments:
        specs.append([{"secondary_y": True}])

    fig = make_subplots(rows=len(treatments), cols=1,
                    subplot_titles = species,
                    specs=specs)
    weight_df = pwg.plastid_weight_sums.reset_index()
    weight_df = weight_df[weight_df['tissue']=='Leaf']
    print(weight_df[weight_df['treatment']=='FeLim'])
    
    print("FeLim ", weight_df[weight_df['treatment']=='FeLim']['TotalPlastidMass'].mean())

    print("Conrol ", weight_df[weight_df['treatment']=='FeLim']['TotalPlastidMass'].mean())

    # print(abc)
    index = 0

    for trmt in treatments:
        index += 1
        score_temp = counts_df[counts_df['treatment']==trmt]
        weight_trmt = weight_df[weight_df['treatment']==trmt]
        weight_ctrl = weight_df[weight_df['treatment'] == 'Control']
        weight_ctrl = weight_ctrl[weight_ctrl['time_stamp']!='00h']
        if verbose: 
            print("Total Plastid Mass Min: ", weight_trmt['TotalPlastidMass'].min())
            print("Total Plastid Mass Max: ", weight_trmt['TotalPlastidMass'].max())

        # Add all the bar traces
        legend = (index==1)
        for var, var_name in method_map.items():#counts_df['variable'].unique():
            fig.add_trace(go.Bar(
                            x=score_temp[score_temp['variable'] == var]['time_stamp'],
                            y=score_temp[score_temp['variable'] == var]['counts'],
                            name=var_name,
                            marker_color=color_map[var],
                            showlegend=legend,
                        ),
                        row=index, col=1,
                        secondary_y=False
            )

        # Add the total plastid mass (138, 147, 186) 113, 121, 153
        fig.add_trace(
            go.Scatter(x=weight_trmt['time_stamp'], y=weight_trmt['TotalPlastidMass'],
                            name=f"MPP {trmt}",
                            # color = 'Treatment',
                            marker_color="rgba(113, 121, 153, 1)",
                            showlegend=legend),
            secondary_y=True,
            row=index, col=1,
        )
        fig.add_trace(
            go.Scatter(x=weight_ctrl['time_stamp'], y=weight_ctrl['TotalPlastidMass'],
                            name="Control",
                            # color = 'Treatment',
                            marker_color="rgba(140, 121, 153, 1)",
                            showlegend=legend,
                            line=dict(dash='dash')),
            secondary_y=True,
            row=index, col=1,
        )
        
    # all treatments 
    if len(treatments) > 1:
        ht, wt = 800, 600
    else:
        ht, wt = 290, 500
        
    fig.update_layout(barmode='group',
                      xaxis_tickangle=-45,
                      height=ht,
                      width=wt,
                      font=dict(
                          family="Arial",
                          size=12
                      ),
                      title = f"{spc.name}",
                      # legend_title_font=dict(size=20),  # Adjust the size of the legend title
                      legend_font=dict(size=12) 
                    )
    if spc.name == 'Sorghum':
        fig.update_layout(legend=dict(
            orientation="h",
            yanchor="top",
            # y=1.4,
            y=-0.25,
            xanchor="left",
            x=0.01,
            # entrywidth=0.2, # change it to 0.3
            # entrywidthmode='fraction',        
            # # itemsizing='constant',  # Use 'constant' to control item size
            # # itemwidth=30,  # Adjust the width of the legend items
            indentation = 0,
            # itemheight=30 
            # sizing="fill",
            # opacity=.7
        ))
    else: 
        fig.update_layout(showlegend=False)


    fig.update_annotations(font_size=12)

    max_plastid = 4.5e-15
    min_plastid = 1.45e-15
    min_count = 4 
    max_count = 162

    # fig.update_yaxes(range=[weight_df['TotalPlastidMass'].min(), weight_df['TotalPlastidMass'].max()], secondary_y=True)
    # fig.update_yaxes(range=[counts_df['counts'].min(), counts_df['counts'].max()], secondary_y=False)
    fig.update_yaxes(range=[min_plastid, max_plastid], secondary_y=True)
    fig.update_yaxes(range=[min_count, max_count], secondary_y=False)
    
    fig.show()
    # print(abc)
    plot_path = os.path.join(csp.results_folder, f"{spc.name}_TotalPlastidMass_RESComps.png")
    pio.write_image(fig, plot_path, scale=6, width=wt, height=ht)

    return 

    print("min plastid: ", weight_df['TotalPlastidMass'].min(), "max plastid: ", weight_df['TotalPlastidMass'].max())
    print("min count: ", counts_df['counts'].min(), "max count: ", counts_df['counts'].max())
    score_only = False
    if score_only:
        fig = px.bar(counts_df, x="time_stamp",
                        y = 'counts',
                        color = 'variable',
                        color_discrete_map=color_map,
                        # facet_col='treatment',
                        facet_row='treatment'
                        # , labels={"rxn_score_I_dist_"+trmt: "Score Dist"}
                        # , category_orders=category_orders
                        , height=800, width=600
                        , facet_row_spacing=0.1
                        , facet_col_spacing=0.03
                        # , height=1000, width=3000
                        )
        fig.update_layout(barmode='group')
        fig.show()
        # print(abc)
        fig = px.histogram(tmm_old_scores, x="value",
                            color = 'variable',
                            facet_col='treatment',
                            facet_row='time_stamp',
                            color_discrete_map=color_map,
                            hover_data=["rxn_ID"]
                            # , labels={"rxn_score_I_dist_"+trmt: "Score Dist"}
                            # , category_orders=category_orders
                            , title = f"{spc.name}"
                            , height=800, width=800
                            , facet_row_spacing=0.1
                            , facet_col_spacing=0.03
                            , marginal="rug"
                            # , height=1000, width=3000
                        )
        fig.update_layout(barmode='group')
        fig.show()

        fig = px.scatter(
            tmm_old_scores[tmm_old_scores['value']>=80], #variable  value
            x='rxn_ID',
            y='value',
            color='variable',
            title='scores: relative sum and max',
            color_discrete_map=color_map,
            # symbol= 'subsystems',
            # color_continuous_scale='icefire',
            # range_color=[-1*mx, mx],
            facet_col='time_stamp',
            facet_row='treatment',
            hover_data=["rxn_ID"]
            # , labels={"rxn_score_I_dist_"+trmt: "Score Dist"}
            # , category_orders=category_orders
            , height=750, width=2500
            , facet_row_spacing=0.1
            , facet_col_spacing=0.03
            # , height=1000, width=3000
        )
        fig.show()

    ## ########################### ########################### ###########################


## ## ----------- OTHER PROCESSING (not used in pipeline) -----------
def compareProteins(aProt, oProt, otherSpc):
    # Check number of modelReactionProteins
    str = ""
    if len(aProt) > len(oProt):
        str += f"\t{otherSpc} has {len(aProt)-len(oProt)} less proteins \n"
    elif len(aProt) < len(oProt):
        str += f"\t{otherSpc} has {len(oProt)-len(aProt)} more proteins\n"
    else:
        for i in range(0, len(aProt)):
            if len(aProt[i]) != len(oProt[i]):
                str += f"\tProt {i+1}: aProt {len(aProt[i])} vs. {len(oProt[i])}\n"
            elif (0 in oProt[i]) and (0 not in aProt[i]):
                str += f"\tProt {i+1}: aProt {aProt[i]} vs. {oProt[i]}\n"

    return str

def checkProteinLists(athalianaModelJson, otherSpcModelJson):
    athalianaModel = dict()
    with open(athalianaModelJson, 'r') as f:
        athalianaModel = json.load(f)
    otherModel = dict()
    with open(otherSpcModelJson, 'r') as f:
        otherModel = json.load(f)
    if not (otherModel and athalianaModel):
        return

    otherSpc = otherModel["name"]

    athaliana_r2p = dict()
    for reaction in athalianaModel["modelreactions"]:
        prot_list = [[len(mrps["feature_refs"])
                        for mrps in mrp["modelReactionProteinSubunits"]]
                            for mrp in reaction["modelReactionProteins"]]
        athaliana_r2p[ reaction["id"]] = prot_list


    otherSpc_r2p = dict()
    for reaction in otherModel["modelreactions"]:
        prot_list = [[len(mrps["feature_refs"])
                        for mrps in mrp["modelReactionProteinSubunits"]]
                            for mrp in reaction["modelReactionProteins"]]
        otherSpc_r2p[ reaction["id"]] = prot_list

    rxns = set()
    for a_rxnID, aProt in athaliana_r2p.items():
        # if len(aProt) >= 3:
        #     print(f"------------------------ {a_rxnID} {len(aProt)}")

        try:
            str = compareProteins(aProt, otherSpc_r2p[a_rxnID], otherSpc)
            if str != "":
                # print(f"{a_rxnID}-{len(aProt)}:\n{str}")
                rxns.add(a_rxnID)
        except KeyError:
            print(f"{a_rxnID} not found in {otherSpc}.")
            c = 0

    return rxns

def rxn2genes(otherSpcModelJson):
    print("get subsystems")
    rxn_subsystems = integrateRoles()
    print("generating file")
    otherModel = dict()
    with open(otherSpcModelJson, 'r') as f:
        otherModel = json.load(f)
    if not otherModel:
        return

    otherSpc = otherModel["name"]

    with open(otherSpc+'_reaction_paralogs.tsv','w') as fh:
        fh.write("reaction_id\tparalog\tenzyme_number\tsubunit_number\tsubsystems\n")
        for reaction in otherModel["modelreactions"]:
            try:
                subsystems = rxn_subsystems[reaction['id'].split('_')[0]]
            except KeyError:
                subsystems = 'Unknown'

            p = 1
            for mrp in reaction["modelReactionProteins"]:
                su = 1
                for mrps in mrp["modelReactionProteinSubunits"]:
                    if len(mrps["feature_refs"]) == 0:
                        fh.write(f"{reaction['id']}\tmissing\t{p}\t{su}\t{subsystems}\n")
                    else:
                        for ftr in mrps["feature_refs"]:
                            fh.write(f"{reaction['id']}\t{ftr.split('/')[-1]}\t{p}\t{su}\t{subsystems}\n")
                    su += 1
                p += 1

def integrateRoles():
    rxn_subsystems = dict()
    print("PS Role Integration ...")
    data = json.load(urlopen(PS_url+PS_tag+PS_json))

    result_list = list()
    subsys_set = set()
    role_set = set()
    for item in data:
        role       = item["role"]
        subsystems = set(item["subsystems"])
        # Get all reactions associated with the role
        reactions = list()
        for r_id in item["reactions"]:
            if r_id in rxn_subsystems:
                rxn_subsystems[r_id] = rxn_subsystems[r_id].union(subsystems)
            else:
                rxn_subsystems[r_id]= subsystems
    return rxn_subsystems

## ## ---------------------------------------------------------------

if __name__ == '__main__':
    #Parameters_test_a, Parameters_QPSI, Parameters_secMeta_TSU
    project = Parameters_QPSI
    generate_reactionScores(project)

