import warnings
warnings.simplefilter(action='ignore', category=Warning)

import plotly.graph_objects as go
import pandas as pa
import numpy as np
from scipy.stats import zscore
import os
from cobra.io import read_sbml_model
import re

import json
from urllib.request import urlopen

from scipy import stats

import plotly.express as px
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import plotly.io as pio
import seaborn as sns

pio.templates.default = "plotly_white" #"none"

from sklearn.decomposition import PCA

import sys
from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from parameters import Parameters

# Identity line
A, B = 1, -1
days = ['02d', '04d', '07d', '14d', '21d']

class bc:
    PROG = '\033[90m'
    RESULT = '\033[92m' #32
    PROMPT = '\033[95m'
    SUBRESULT = '\033[34m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'
    #96

class ResultProcessingHelper:
    def __init__(self, parameters, spc='Sorghum', project="QPSI"):
        self.spc = spc
        self.project = parameters.project
        self.control_name = parameters.control_id
        self.group_cols = parameters.group_columns
        self.trmt_col = parameters.trmt_colmn
        self.value = parameters.value_column

        # Write the results to this folder
        self.results_folder = parameters.results_folder

        self.rxn_scores_file = parameters.rxn_scores_paths[spc]
        self.RNASeq_file_path = parameters.rnaseq_paths[spc]

        self.output_folder = os.path.join(parameters.results_folder,"..","processing_results")
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        folder = parameters.json_files_folder
        if self.project == "QPSI":
            spc_model = "sbicolor_3.1.1_Thylakoid_Reconstruction_ComplexFix_070224.xml" if spc == 'Sorghum' else "ptrich_4.1_Thylakoid_Reconstruction_ComplexFix_070224.xml"

            if 'plastid' in param.json_files_folder:
                spc_model = "sbicolor_3.1.1_plastid_Thylakoid_Reconstruction_ComplexFix_070224_noADP.xml" if spc == 'Sorghum' else "ptrich_4.1_plastid_Thylakoid_Reconstruction_ComplexFix_070224_noADP.xml"
            self.sbmlModel = f"{folder}/{spc_model}"

        
        self.cols_values = None
        # self.quantiles = [i/10 for i in range(0, 11)]
        self.quantiles = [i*5/100 for i in range(0, 21)]
        self.labels = [int(i*100) for i in self.quantiles[:-1]]

        self.rxn_scores_df = None
        self.rna_df = None
        self.all_features_dict = None
        self.fluxes_df = None
        self.rxn_ftrExpression_df = None
        self.significant_reactions_features = None

        self.score_perct = 95

        self.sig_subsystems = pa.DataFrame()

## ----------------------- Helper methods ------------------------------------------
def apply_literal_eval(row, col='subsystems'):
    if isinstance(row[col], str):
        import ast
        if "Calvin-Benson-Bassham_cycle_in_plants" in row[col]:
            return ["Calvin-Benson-Bassham_cycle_in_plants"]
        else:
            return ast.literal_eval(row[col])
    else:
        return []

## ---------------------------------------------------------------------------------

## ----------------------- PLOT REACTION SCORE SCATTER -----------------------------
def edit_rxn_scores_df(RPHelper, rxn_scores_file=''):
    print(" --------- Loading reactions score DF --------- ")
    if 'day' in RPHelper.group_cols:
            RPHelper.group_cols[RPHelper.group_cols.index('day')] = 'time_stamp'

    if not rxn_scores_file: 
        rxn_scores_file = RPHelper.rxn_scores_file

    sep = '\t' if '.tsv' in rxn_scores_file else ','

    df = pa.read_csv(rxn_scores_file, sep=sep)
    print(df.head())

    if ('norm_value' not in df.columns):
        gcs = RPHelper.group_cols.copy()
        gcs.remove(RPHelper.trmt_col)

        value = 'reaction_score'
        df[value] = df[value].astype('float')
        df[value] = np.log(df[value])
        # Group the scores DF
        groups = df.groupby(gcs)
        # computes group-wise mean/std,
        # then auto broadcasts to size of group chunk
        min = groups[value].transform("min")
        max = groups[value].transform("max")
        df['norm_'+value] = (df[value] - min) / (max - min)

    # Remove unnecessary columns
    keep = {'reaction_id', 'flexibility', 'bind', 'isTrans', 'subsystems', 'features', 'Gene_ID', 'flexibility', 'value', 'norm_value', 'rxn_score_I_dist', 'rxn_dist_quantile', 'spc'}
    keep.update(RPHelper.group_cols)
    df.drop(set(df.columns)-keep, inplace=True, axis=1)

    df['value'] = df['norm_value']
    # Rename time stamp column and keep only used time points
    if RPHelper.project == 'QPSI':
        print(df.columns)
        df.rename(columns={'time_stamp':'day'}, inplace=True)
        if 'time_stamp' in RPHelper.group_cols:
            RPHelper.group_cols[RPHelper.group_cols.index('time_stamp')] = 'day'
        
        df = df[df['day'].isin(days)]


    RPHelper.cols_values = dict()
    for col in RPHelper.group_cols:
        RPHelper.cols_values[col] = sorted(list(df[col].unique()))

    try:
        RPHelper.cols_values[RPHelper.trmt_col].remove(RPHelper.control_name)
    except KeyError:
        print("Group columns did not load properly")

    return df


def sig_subsystems_table(RPHelper, distance_df, dist):
    sig = distance_df[(np.abs(distance_df['rxn_score_I_dist']) >= dist)]
    sig = sig[RPHelper.group_cols+['subsystems', 'class']]
    sig = sig.drop_duplicates()

    ################ remove
    if RPHelper.project != "QPSI":
        tissues = ["top_leaf", "internode_leaf_5", "internode_5", "roots"]#["top_leaf", "internode_leaf_5", "internode_5", "roots"]
        trmts = ["N20", "N40", "N60", "N80"]#["N80", "N40", "N20"]
        days = ["vegetative", "boot", "anthesis", "grain_fill"]#["anthesis", "vegetative"]

        sig = sig[(sig['tissue'].isin(tissues)) &
                  (sig['treatment'].isin(trmts)) &
                  (sig['growth_stage'].isin(days))]

    sig.to_csv(f'{RPHelper.output_folder}{RPHelper.spc}_sig_subsystems.csv', index=False)
    ################

    if RPHelper.project == 'QPSI':
        sig = sig[['tissue', 'treatment', 'day', 'subsystems', 'class']]
        sig['day'] = sig['day'].transform(lambda ts:ts if ts in ['14d', '21d'] else '0'+ts)
        sig.sort_values(by=['tissue', 'treatment', 'day'], inplace=True)
    else:
        sig.sort_values(by=['tissue', 'growth_stage', 'treatment'], inplace=True)
    # print(sig.T)
    sig['in'] = 'X'
    val = ['in']#,
    ind = ['subsystems', 'class']
    #    change from long to wide
    if RPHelper.project == 'QPSI':
        cols = ['tissue', 'treatment', 'day']
    else:
        cols = ['tissue', 'growth_stage', 'treatment']
    sig_piv = sig.pivot(index=ind, columns=cols, values=val)

    # print(sig_piv)
    sig_piv.to_csv(f'{RPHelper.output_folder}{RPHelper.spc}_sig_subsystems_in.csv')

def plot_rxn_scores_subsystem(RPHelper, cmprt='', verbose=False, heatmap=False):

    if verbose: print(RPHelper.rxn_scores_df.columns)
    
    rxn_scores_df = RPHelper.rxn_scores_df[['rxn_ID', 'bind', 'norm_value', \
        'rxn_score_I_dist', 'subsystems', 'rxn_dist_quantile', 'features']\
        +RPHelper.group_cols].copy()
    
    rxn_scores_df.rename(columns={'norm_value':'value'}, inplace=True)


    ctrl = RPHelper.control_name
    cols_indices = [i for i in range(0, len(RPHelper.group_cols))]
    if RPHelper.project == 'QPSI': cols_indices.reverse()
    
    cols_indices.remove(RPHelper.group_cols.index(RPHelper.trmt_col))
    percentile = RPHelper.score_perct

    # "percentile" distance to the identity line across all treatments, time points  
    dist = np.abs(rxn_scores_df['rxn_score_I_dist'])\
        .describe([percentile/100])[str(percentile)+'%']
    # group"percentile" distance to the identity line across all treatments, time points
    dist_dict = dict()
    group_cols = [RPHelper.group_cols[idx] for idx in cols_indices]
    groups = rxn_scores_df.groupby(group_cols)
    for name, group in groups:
        name = (name[0], name[0]) if name[-1] == '' else (name[0], name[-1])
        dist_dict[name] = np.abs(group['rxn_score_I_dist'])\
            .describe([percentile/100])[str(percentile)+'%']
        # dist_dict[name] = dist

    if verbose: print(f"Distance at {percentile} percentile: {dist}")
    



    rxn_scores_df['subsystems'] = rxn_scores_df.apply(lambda row: apply_literal_eval(row), axis=1)
    rxn_scores_df = rxn_scores_df.explode('subsystems')
    sys_class = get_subsysClass(list(rxn_scores_df['subsystems'].unique()))


    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants_and_prokaryotes','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_',' ')
    # # if you need to remove some subsystems
    # remove_sub = [...]
    # rxn_scores_df = rxn_scores_df[~rxn_scores_df['subsystems'].isin(remove_sub)]
    rxn_scores_df['class'] = rxn_scores_df['subsystems'].map(sys_class)

    # write subsytems to the file
    temp_df = rxn_scores_df[RPHelper.group_cols+['subsystems']]
    temp_df['source'] = 'rxnScore'
    RPHelper.sig_subsystems = pa.concat(
        [RPHelper.sig_subsystems, temp_df],
        ignore_index=True)


    col = [RPHelper.trmt_col]
    val = ['value', 'rxn_score_I_dist', 'rxn_dist_quantile']#,
    ind = set(rxn_scores_df.columns) - set().union(col, val)
    #    change from long to wide
    distance_df = rxn_scores_df.pivot(index=ind, columns=col, values=val)
    distance_df.columns = distance_df.columns.get_level_values(0) + '_' +  distance_df.columns.get_level_values(1)
    distance_df = distance_df.reset_index()


    if verbose: print(distance_df.columns)
    treatments = RPHelper.cols_values[RPHelper.trmt_col].copy()

    wt, ht = 950, 450

    names = {'value_'+trmt : trmt for trmt in treatments}
    names[f'value_{RPHelper.control_name}'] = RPHelper.control_name
    distance_df.rename(columns=names, inplace=True)
    
    if 'time_stamp' in  RPHelper.cols_values:
        distance_df.sort_values(by=['time_stamp'], inplace=True)
    else:
        distance_df.sort_values(by=['day'], inplace=True)


    for trmt in ['FeLim']: # treatments
        # Normalized scores scatter plots
        mx = distance_df["rxn_score_I_dist_"+trmt].max()
        mx = max(np.abs(distance_df["rxn_score_I_dist_"+trmt].min()), mx)
        if verbose: print(f"--> {trmt} max: {mx} {cmprt}")
        
        mx = 0.25
        title = f"{RPHelper.spc} {trmt} -- Normalized reaction scores."

        category_orders = {} if RPHelper.project == "secMeta" \
                            else {"day": days,
                                    "tissue": RPHelper.cols_values['tissue'].copy()} # tissues}

        fig = px.scatter(
                distance_df,
                x=ctrl,
                y=trmt,
                color="rxn_score_I_dist_"+trmt,
                title=title,
                # symbol= 'subsystems',
                # size = 'size',
                color_continuous_scale='icefire',
                range_color=[-1*mx, mx],
                facet_col=RPHelper.group_cols[cols_indices[0]],
                facet_row=RPHelper.group_cols[cols_indices[-1]],
                hover_data=["rxn_ID", "bind"],
                labels={"rxn_score_I_dist_"+trmt: "Score Dist"},
                category_orders=category_orders
                , height=ht, width=wt
                , facet_row_spacing=0.08
                , facet_col_spacing=0.03
            )
        reference_line_x_range = np.array([0, 1])
        fig.update_traces(marker=dict(size=6))


        col_values = RPHelper.cols_values[RPHelper.group_cols[cols_indices[0]]]
        row_values = RPHelper.cols_values[RPHelper.group_cols[cols_indices[-1]]]

        num_unique = 0
        for row_idx, row_figs in enumerate(fig._grid_ref):
            for col_idx, col_fig in enumerate(row_figs):
                # Add the greater variability subsystems reactions
                col_val = col_values[col_idx] if isinstance(col_values, list) else col_values
                row_val = row_values[(row_idx+1)*-1] if isinstance(row_values, list) else row_values

                trmt_df = distance_df[(np.abs(distance_df['rxn_score_I_dist_'+trmt]) >= \
                                        dist_dict[(col_val, row_val)])
                                      & (distance_df[RPHelper.group_cols[cols_indices[0]]] == col_val)
                                      & (distance_df[RPHelper.group_cols[cols_indices[-1]]] == row_val)]

                fig = fig.add_trace(go.Scatter(x=trmt_df[ctrl],
                                        y=trmt_df[trmt],
                                        # x=sys_trmt_1[trmt],
                                        # y=sys_trmt_1[ctrl],
                                        mode='markers',
                                        marker=dict(symbol='circle', size=9, color='black'),
                                        showlegend=False),
                                        row=row_idx+1, col=col_idx+1)


        # Add identity line to the plots
        fig =fig.add_trace(go.Scatter(x=reference_line_x_range,
                                        y=A*reference_line_x_range,
                                        marker_color='rgba(65, 65, 65, .4)',
                                        showlegend=False),
                                        row='all', col='all')

        fig.update_xaxes(range=[-0.01,0.7], showticklabels=True)
        fig.update_yaxes(range=[-0.01,0.7], showticklabels=True)
        
        fig.update_layout(
            showlegend=False,
            font=dict(
                family="arial",#Courier New, Arial
                size=14,
            ),
            legend=dict(
                font=dict(
                    family="arial",
                    size=14,
                    color="black"
                )
            ),
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))


        fig.data = [fig.data[i] for i in reversed(range(len(fig.data)))]

        plot_path = RPHelper.results_folder+f"scatter_plots_{RPHelper.spc}_{trmt}.png"
        pio.write_image(fig, plot_path, scale=2, width=wt, height=ht)

    compute_var_std(rxn_scores_df, RPHelper.trmt_col, RPHelper.spc)

def compute_var_std(input_df, trmt_col, spc='', cprmt='', var_column='rxn_score_I_dist', verbose=False, cols_values='', group_cols=''):
    print(" --------- Processing RIPTiDe results * Fluxes scatter std line plots --------- ")

    if not os.path.exists(RPHelper.results_folder+'std_figures/'):
        os.makedirs(RPHelper.results_folder+'std_figures/')

    print(plt.style.available)
    plt.style.use('seaborn-v0_8-bright')
    plt.rcParams['font.family'] = "Times New Roman"
    # plt.rcParams['font.serif'] = 'Ubuntu'
    # plt.rcParams['font.monospace'] = 'Ubuntu Mono'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 10
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.titlesize'] = 12


    # sort_dict = {'growth_stage':["vegetative", "boot", "anthesis", "grain_fill"],
    #              'treatment': ['N80', 'N60', 'N40', 'N20'],
    #              'tissue': ["top_leaf", "internode_leaf_5", "internode_5", "stem_top", "roots"]}
    # sort_dict = {'growth_stage':["vegetative", "anthesis"],
    #              'treatment': ['N80', 'N60', 'N40', 'N20'],
    #              'tissue': ["top_leaf", "internode_leaf_5", "internode_5", "roots"]}
    if not cols_values:
        sort_dict = RPHelper.cols_values
    else: 
        sort_dict = cols_values
        if trmt_col not in sort_dict:
             sort_dict[trmt_col] = RPHelper.cols_values[trmt_col]

    if not group_cols:
        group_cols = RPHelper.group_cols #['tissue', 'growth_stage', 'treatment']
    # else: 


    if 'tissue' in input_df:
        df = input_df[input_df['tissue']=='Leaf'].copy()
    else: 
        df = input_df.copy()
    if cprmt != '':
        print(f"cprmt {df.shape}")
        df = df[(df['rxn_ID'].str.contains(cprmt))]
        print(f"cprmt {df.shape}")

    # Remove the tissue constraint
    df = df[(df[trmt_col].isin(sort_dict[trmt_col]))]
    print(df.head())
    # ----- Compute var and std for plots
    if verbose: print(f"upregulated std")
    up_std = df[df[var_column] >= 0][group_cols+[var_column]].\
        groupby(group_cols).var().unstack(1)
    up_std.columns = up_std.columns.get_level_values(1)
    print("-----------> up")
    print(up_std.head())
    print(sort_dict)
    print(group_cols)

    if not cols_values:
        up_trans = up_std.T
        print(up_trans)
        up_trans.columns = up_trans.columns.get_level_values(1)
    else: 
        up_trans = up_std

    up_trans = up_trans[sort_dict[group_cols[-1]]]
    if not cols_values:
        up_trans = up_trans.T
    print(up_trans.head())
    # print(abc)


    if verbose: print(f"downregulated std")
    down_std = df[df[var_column] <= 0][group_cols+[var_column]].\
        groupby(group_cols).var().unstack(1)
    down_std.columns = down_std.columns.get_level_values(1)
    
    if not cols_values:
        down_trans = down_std.T
        down_trans.columns = down_trans.columns.get_level_values(1)
    else: 
        down_trans = down_std
    down_trans = down_trans[sort_dict[group_cols[-1]]]
    
    if not cols_values:
        down_trans = down_trans.T
    print("-----------> up")
    print(down_trans.head())
    # print(abc)

    upmax = up_trans.select_dtypes(include=[np.number]).max().max()
    downmax = down_trans.select_dtypes(include=[np.number]).max().max()
    ymax = max(upmax, downmax)
    if verbose: print(f"---------------  Maximum STD {ymax}")
    print(f"---------------  Maximum STD {ymax}")
    ymax = 0.000244
    if cols_values: ymax = 0.0008

    print(up_trans.head())
    print(group_cols)

    res_label= {"relab": r"$r_{s}^{r} \; here$",
                "obj": r"$r_{s} \; here$"}
    res_labels = {r+"_"+s: res_label[r].replace('here', s) for r in ['obj', 'relab']
                        for s in ['Poplar', 'Sorghum']} 
    for trmt in sort_dict[group_cols[1]]:
        col = f'{trmt}' #if not cols_values else res_labels[trmt]
        if verbose: print(f"Generating figure for {col}")
        width, height = 3.5, 2
        fig, ax = plt.subplots(figsize=(width, height))
        ax.plot(up_trans[col], label='Up-regulated', color='chocolate')
        ax.plot(down_trans[col], label='Down-regulated', color='steelblue')
        # ax.legend(fancybox=True)
        # leg=plt.legend(loc='best', numpoints=1, fancybox=True)
        ax.grid()
        ax.ticklabel_format(style='sci', axis='y', useOffset=True)
        # plt.ticklabel_format(style='sci', useOffset=True)
        plt.xlabel('Time point')
        plt.ylabel(r'$\sigma$')
        ax.set_ylim([0, ymax+0.00005])
       
        
        plt.legend(frameon=True, loc='best', facecolor='white', framealpha=.8, fancybox=True)
        
        title = f"{col.replace('_', ' ')} {cprmt}" if not cols_values else res_labels[trmt]
        fig.suptitle(title)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.savefig(RPHelper.results_folder+f'std_figures/{spc}_{col}_Leaf.png', bbox_inches='tight', dpi=400)

## ----------------------- PAPER FIGURES - PROJECT SPECIFIC -----------------------
def get_subsysClass(subsystems, verbose=False):
    pathways_class_dict = dict()
    PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
    PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
    PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"
    PS_json_data = json.load(urlopen(PS_url+PS_tag+PS_json))
    if verbose: print("LOADED ROLES")
    notProcessed = subsystems.copy()

    for item in PS_json_data:
        all_classes = item["classes"].keys()
        for cls in all_classes:
            for pathway in item["classes"][cls]:
                if pathway in notProcessed:
                    new_pathway = pathway.replace('_in_plants','')
                    new_pathway = new_pathway.replace('_',' ')
                    pathways_class_dict[new_pathway] = cls
                    pathways_class_dict[f"Z{cls}"] = cls
                    notProcessed.remove(pathway)
                if not notProcessed: break

    pathways_class_dict['Ferredoxins'] = 'Other'
    for fr in ["1", "2", "3", "C1", "C2"]:
        pathways_class_dict[f'Ferredoxin {fr}'] = 'Other'
    return pathways_class_dict


def all_rxn_scores(RPHelper, rxn_scores_df, cmprt='', verbose=False, heatmap=False):
    rxn_scores_df.rename(columns={'norm_value':'value'}, inplace=True)
    group_cols = ['spc', 'res']
    group_cols = ['day', 'res_spc'] #spc_res
    res_spc = [r+'_'+s for r in ['obj', 'relab']
                        for s in ['Poplar', 'Sorghum']] 

    spc_res = [s+'_'+r for s in ['Poplar', 'Sorghum']
                        for r in ['relab', 'obj']] 
    cols_values = {'spc': ['Poplar', 'Sorghum'], 
                    'res': ['relab', 'obj'], 
                    'spc_res': spc_res, #['Poplar_relab', 'Poplar_obj', 'Sorghum_relab', 'Sorghum_obj'],
                    'res_spc' : res_spc, 
                    'day': RPHelper.cols_values['day']}

    res_label= {"relab": r"$r_{s}^{r} \; here$",
                "obj": r"$r_{s} \; here$"}
    res_labels = {r+"_"+s: res_label[r].replace('here', s) for r in ['obj', 'relab']
                        for s in ['Poplar', 'Sorghum']} 
    print(res_labels)




    ctrl = RPHelper.control_name
    cols_indices = [i for i in range(0, len(group_cols))]
    
    percentile = RPHelper.score_perct

    dist = np.abs(rxn_scores_df['rxn_score_I_dist'])\
        .describe([percentile/100])[str(percentile)+'%']

    dist_dict = dict()

    group_cols = [group_cols[idx] for idx in cols_indices]
    groups = rxn_scores_df.groupby(group_cols)


    for name, group in groups:
        print(name)
        name = (name[0], name[0]) if name[-1] == '' else (name[0], name[-1])
        dist_dict[name] = np.abs(group['rxn_score_I_dist'])\
            .describe([percentile/100])[str(percentile)+'%']
        # dist_dict[name] = dist
    print(dist_dict)


    # print(abc)
    if verbose: print(f"Distance at {percentile} percentile: {dist}")
    print(f"Distance at {percentile} percentile: {dist}")



    rxn_scores_df['subsystems'] = rxn_scores_df.apply(lambda row: apply_literal_eval(row), axis=1)
    rxn_scores_df = rxn_scores_df.explode('subsystems')
    sys_class = get_subsysClass(list(rxn_scores_df['subsystems'].unique()))

    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants_and_prokaryotes','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_',' ')
    rxn_scores_df['class'] = rxn_scores_df['subsystems'].map(sys_class)



    col = [RPHelper.trmt_col]
    val = ['value', 'rxn_score_I_dist', 'rxn_dist_quantile']#,
    ind = set(rxn_scores_df.columns) - set().union(col, val)
    #    change from long to wide
    distance_df = rxn_scores_df.pivot(index=ind, columns=col, values=val)
    distance_df.columns = distance_df.columns.get_level_values(0) + '_' +  distance_df.columns.get_level_values(1)
    distance_df = distance_df.reset_index()
    print(distance_df.shape)



    if verbose: print(distance_df.columns)

    treatments = RPHelper.cols_values[RPHelper.trmt_col].copy()
    # print(f"{treatments} {ctrl}")
    # treatments.remove(ctrl)
    if (RPHelper.project != 'QPSI') and (cmprt != 'all'): tmts = ['N20']
    else: tmts = treatments


    if RPHelper.project == 'hAlpha':
        wt, ht = 1400, 1000
    elif RPHelper.project == 'QPSI':
        wt, ht = 1400, 600
        wt, ht = 800, 600 # the best
        wt, ht = 1000, 800 # the best
    else:
        wt, ht = 1650, 600

    names = {'value_'+trmt : trmt for trmt in tmts}
    names[f'value_{RPHelper.control_name}'] = RPHelper.control_name
    distance_df.rename(columns=names, inplace=True)
    if 'time_stamp' in  RPHelper.cols_values:
        distance_df.sort_values(by=['time_stamp'], inplace=True)
    else:
        distance_df.sort_values(by=['day'], inplace=True)

    # distance_df = distance_df[distance_df['tissue']=='Leaf']
    # distance_df['size'] = 0.1
    for trmt in ['FeLim']: # tmts:#RPHelper.cols_values[RPHelper.trmt_col]:
    # for trmt in ['FeLim', 'ZnEx']:#RPHelper.cols_values[RPHelper.trmt_col]:
        # Normalized scores scatter plots
        mx = distance_df["rxn_score_I_dist_"+trmt].max()
        mx = max(np.abs(distance_df["rxn_score_I_dist_"+trmt].min()), mx)
        if verbose: print(f"--> {trmt} max: {mx} {cmprt}")
        print(f"--> {trmt} max: {mx} {cmprt}")
        mx = 0.25
        title = f"{RPHelper.spc} {trmt} -- Normalized reaction scores."

        category_orders = {} if RPHelper.project == "secMeta" \
                            else {"day": days, 
                                    "tissue": RPHelper.cols_values['tissue'].copy(), 
                                    "res_spc": res_spc} # tissues}

        lebels = res_labels
        lebels["rxn_score_I_dist_"+trmt] =  "Score Dist"
        print(lebels)
        # print(abc)
        fig = px.scatter(
                distance_df,
                # x=trmt,
                # y=ctrl,
                x=ctrl,
                y=trmt,
                color="rxn_score_I_dist_"+trmt,
                title=title,
                # symbol= 'subsystems',
                # size = 'size',
                color_continuous_scale='icefire',
                range_color=[-1*mx, mx],
                facet_col=group_cols[cols_indices[0]],
                facet_row=group_cols[cols_indices[-1]],
                hover_data=["rxn_ID", "bind"],
                labels=lebels,
                category_orders=category_orders
                , height=ht, width=wt
                , facet_row_spacing=0.02
                , facet_col_spacing=0.02
                # , height=1000, width=3000
            )
        reference_line_x_range = np.array([0, 1])
        fig.update_traces(marker=dict(size=6))


        col_values = cols_values[group_cols[cols_indices[0]]]
        row_values = cols_values[group_cols[cols_indices[-1]]]

        num_unique = 0
        for row_idx, row_figs in enumerate(fig._grid_ref):
            for col_idx, col_fig in enumerate(row_figs):
                # Add the greater variability subsystems reactions
                col_val = col_values[col_idx] if isinstance(col_values, list) else col_values
                row_val = row_values[(row_idx+1)*-1] if isinstance(row_values, list) else row_values

                trmt_df = distance_df[(np.abs(distance_df['rxn_score_I_dist_'+trmt]) >= \
                                        dist_dict[(col_val, row_val)])
                                      & (distance_df[group_cols[cols_indices[0]]] == col_val)
                                      & (distance_df[group_cols[cols_indices[-1]]] == row_val)]


                fig = fig.add_trace(go.Scatter(x=trmt_df[ctrl],
                                        y=trmt_df[trmt],
                                        # x=sys_trmt_1[trmt],
                                        # y=sys_trmt_1[ctrl],
                                        mode='markers',
                                        marker=dict(symbol='circle', size=9, color='black'),
                                        showlegend=False),
                                        row=row_idx+1, col=col_idx+1)

        # Add identity line to the plots
        fig =fig.add_trace(go.Scatter(x=reference_line_x_range,
                                        y=A*reference_line_x_range,
                                        marker_color='rgba(65, 65, 65, .4)',
                                        showlegend=False),
                                        row='all', col='all')

        # fig.update_xaxes(range=[-0.01, 0.62], showticklabels=True)
        # fig.update_yaxes(range=[-0.01, 0.62])#, showticklabels=True)

        # fig.update_xaxes(range=[-0.01,0.5])#, showticklabels=True)
        # fig.update_yaxes(range=[-0.01,0.5])#, showticklabels=True)
        fig.update_xaxes(range=[-0.01,0.65])#, showticklabels=True)
        fig.update_yaxes(range=[-0.01,0.65])#, showticklabels=True)
        fig.update_layout(
            showlegend=False,
            font=dict(
                family="arial",#Courier New, Arial
                size=14,
            ),
            legend=dict(
                font=dict(
                    family="arial",
                    size=14,
                    color="black"
                )
            ),
        )
        fig.for_each_annotation(lambda a: a.update(text=res_labels[a.text.split("=")[1]] \
                                                    if 'res' in a.text.split("=")[0] \
                                                    else a.text.split("=")[1] \
                                    ))


        fig.data = [fig.data[i] for i in reversed(range(len(fig.data)))]
        fig.show()
        plot_path = RPHelper.results_folder+f"scatter_plots_{RPHelper.spc}_{trmt}.png"
        pio.write_image(fig, plot_path, scale=4, width=wt, height=ht)

    compute_var_std(rxn_scores_df, RPHelper.trmt_col, RPHelper.spc , cols_values=cols_values , group_cols=group_cols)



def all_rxn_scores_spc(RPHelper, rxn_scores_df, cmprt='', verbose=False, heatmap=False):
    rxn_scores_df.rename(columns={'norm_value':'value'}, inplace=True)
    rxn_scores_df = rxn_scores_df[rxn_scores_df['treatment'].isin(['Control', 'FeLim'])]


    fig_rows = ['Control', 'FeLim']

    group_cols = ['treatment', 'res']
    group_cols = ['day', 'res_trmt'] #spc_res
    res_trmt = [r+'_'+s for r in ['obj', 'relab']
                        for s in fig_rows] 

    trmt_res = [s+'_'+r for s in fig_rows
                        for r in ['relab', 'obj']] 
    cols_values = {'treatment': fig_rows, 
                    'res': ['relab', 'obj'], 
                    'trmt_res': trmt_res, #['Poplar_relab', 'Poplar_obj', 'Sorghum_relab', 'Sorghum_obj'],
                    'res_trmt' : res_trmt, 
                    'day': RPHelper.cols_values['day']}

    res_label= {"relab": r"$r_{s}^{r} \; here$",
                "obj": r"$r_{s} \; here$"}
    res_labels = {r+"_"+s: res_label[r].replace('here', s) for r in ['obj', 'relab']
                        for s in fig_rows} 
    print(res_labels)




    ctrl = RPHelper.control_name
    cols_indices = [i for i in range(0, len(group_cols))]
    
    percentile = RPHelper.score_perct

    dist = np.abs(rxn_scores_df['rxn_score_I_dist'])\
        .describe([percentile/100])[str(percentile)+'%']

    dist_dict = dict()

    group_cols = [group_cols[idx] for idx in cols_indices]
    groups = rxn_scores_df.groupby(group_cols)


    for name, group in groups:
        print(name)
        name = (name[0], name[0]) if name[-1] == '' else (name[0], name[-1])
        dist_dict[name] = np.abs(group['rxn_score_I_dist'])\
            .describe([percentile/100])[str(percentile)+'%']
        # dist_dict[name] = dist
    print(dist_dict)


    # print(abc)
    if verbose: print(f"Distance at {percentile} percentile: {dist}")
    print(f"Distance at {percentile} percentile: {dist}")



    rxn_scores_df['subsystems'] = rxn_scores_df.apply(lambda row: apply_literal_eval(row), axis=1)
    rxn_scores_df = rxn_scores_df.explode('subsystems')
    sys_class = get_subsysClass(list(rxn_scores_df['subsystems'].unique()))

    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants_and_prokaryotes','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants','')
    rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_',' ')
    rxn_scores_df['class'] = rxn_scores_df['subsystems'].map(sys_class)



    col = [RPHelper.trmt_col]
    val = ['value', 'rxn_score_I_dist', 'rxn_dist_quantile']#,
    ind = set(rxn_scores_df.columns) - set().union(col, val)
    #    change from long to wide
    distance_df = rxn_scores_df.pivot(index=ind, columns=col, values=val)
    distance_df.columns = distance_df.columns.get_level_values(0) + '_' +  distance_df.columns.get_level_values(1)
    distance_df = distance_df.reset_index()
    print(distance_df.shape)



    if verbose: print(distance_df.columns)

    treatments = RPHelper.cols_values[RPHelper.trmt_col].copy()

    if (RPHelper.project != 'QPSI') and (cmprt != 'all'): tmts = ['N20']
    else: tmts = treatments

    if RPHelper.project == 'hAlpha':
        wt, ht = 1400, 1000
    elif RPHelper.project == 'QPSI':
        wt, ht = 1400, 600
        wt, ht = 800, 600 # the best
        wt, ht = 1000, 800 # the best
    else:
        wt, ht = 1650, 600

    names = {'value_'+trmt : trmt for trmt in tmts}
    names[f'value_{RPHelper.control_name}'] = RPHelper.control_name
    distance_df.rename(columns=names, inplace=True)
    if 'time_stamp' in  RPHelper.cols_values:
        distance_df.sort_values(by=['time_stamp'], inplace=True)
    else:
        distance_df.sort_values(by=['day'], inplace=True)

    # distance_df = distance_df[distance_df['tissue']=='Leaf']
    # distance_df['size'] = 0.1
    for trmt in ['FeLim']: # tmts:#RPHelper.cols_values[RPHelper.trmt_col]:
    # for trmt in ['FeLim', 'ZnEx']:#RPHelper.cols_values[RPHelper.trmt_col]:
        # Normalized scores scatter plots
        mx = distance_df["rxn_score_I_dist_"+trmt].max()
        mx = max(np.abs(distance_df["rxn_score_I_dist_"+trmt].min()), mx)
        if verbose: print(f"--> {trmt} max: {mx} {cmprt}")
        print(f"--> {trmt} max: {mx} {cmprt}")
        mx = 0.25
        title = f"{RPHelper.spc} {trmt} -- Normalized reaction scores."

        category_orders = {} if RPHelper.project == "secMeta" \
                            else {"day": days,
                                    "tissue": RPHelper.cols_values['tissue'].copy(), 
                                    "res_spc": res_spc} # tissues}

        lebels = res_labels
        lebels["rxn_score_I_dist_"+trmt] =  "Score Dist"
        print(lebels)
        # print(abc)
        fig = px.scatter(
                distance_df,
                # x=trmt,
                # y=ctrl,
                x=ctrl,
                y=trmt,
                color="rxn_score_I_dist_"+trmt,
                title=title,
                # symbol= 'subsystems',
                # size = 'size',
                color_continuous_scale='icefire',
                range_color=[-1*mx, mx],
                facet_col=group_cols[cols_indices[0]],
                facet_row=group_cols[cols_indices[-1]],
                hover_data=["rxn_ID", "bind"],
                labels=lebels,
                category_orders=category_orders
                , height=ht, width=wt
                , facet_row_spacing=0.02
                , facet_col_spacing=0.02
                # , height=1000, width=3000
            )
        reference_line_x_range = np.array([0, 1])
        fig.update_traces(marker=dict(size=6))


        col_values = cols_values[group_cols[cols_indices[0]]]
        row_values = cols_values[group_cols[cols_indices[-1]]]

        num_unique = 0
        for row_idx, row_figs in enumerate(fig._grid_ref):
            for col_idx, col_fig in enumerate(row_figs):
                # Add the greater variability subsystems reactions
                col_val = col_values[col_idx] if isinstance(col_values, list) else col_values
                row_val = row_values[(row_idx+1)*-1] if isinstance(row_values, list) else row_values
                trmt_df = distance_df[(np.abs(distance_df['rxn_score_I_dist_'+trmt]) >= \
                                        dist_dict[(col_val, row_val)])
                                      & (distance_df[group_cols[cols_indices[0]]] == col_val)
                                      & (distance_df[group_cols[cols_indices[-1]]] == row_val)]


                fig = fig.add_trace(go.Scatter(x=trmt_df[ctrl],
                                        y=trmt_df[trmt],
                                        # x=sys_trmt_1[trmt],
                                        # y=sys_trmt_1[ctrl],
                                        mode='markers',
                                        marker=dict(symbol='circle', size=9, color='black'),
                                        showlegend=False),
                                        row=row_idx+1, col=col_idx+1)

        # Add identity line to the plots
        fig =fig.add_trace(go.Scatter(x=reference_line_x_range,
                                        y=A*reference_line_x_range,
                                        marker_color='rgba(65, 65, 65, .4)',
                                        showlegend=False),
                                        row='all', col='all')

        fig.update_xaxes(range=[-0.01, 0.62], showticklabels=True)
        fig.update_yaxes(range=[-0.01, 0.62])#, showticklabels=True)

        fig.update_layout(
            showlegend=False,
            font=dict(
                family="arial",#Courier New, Arial
                size=14,
            ),
            legend=dict(
                font=dict(
                    family="arial",
                    size=14,
                    color="black"
                )
            ),
        )
        fig.for_each_annotation(lambda a: a.update(text=res_labels[a.text.split("=")[1]] \
                                                    if 'res' in a.text.split("=")[0] \
                                                    else a.text.split("=")[1] \
                                    ))


        fig.data = [fig.data[i] for i in reversed(range(len(fig.data)))]
        fig.show()
        plot_path = RPHelper.results_folder+f"scatter_plots_{RPHelper.spc}_{trmt}.png"
        pio.write_image(fig, plot_path, scale=2, width=wt, height=ht)

    compute_var_std(rxn_scores_df, RPHelper.trmt_col, RPHelper.spc , cols_values=cols_values , group_cols=group_cols)


if __name__ == '__main__':

    params = Parameters()
    spc = "Sorghum"
    RPHelper = ResultProcessingHelper(params, spc)

    # ## ## ## Get reaction scores and the identity line distance
    RPHelper.rxn_scores_df = edit_rxn_scores_df(RPHelper)
    print(RPHelper.rxn_scores_df.head())


    rxn_scores_subsys = True

    if rxn_scores_subsys:
        plot_rxn_scores_subsystem(RPHelper, 'all')


    all_df = pa.DataFrame()
    for spc in ['Poplar', 'Sorghum']: 
        for msr in ['relab', 'obj']: 
            if msr == 'obj':
                rxn_file = RPHelper.results_folder+f"{spc}_objective_abundance_{RPHelper.control_name}.tsv"
            else: 
                rxn_file = RPHelper.results_folder+f"{spc}_relab_rxn_scores_tmm.csv"
            
            scores_df = edit_rxn_scores_df(RPHelper, rxn_file)
            # if 'relab' in RPHelper.rxn_scores_file.lower():
            scores_df = scores_df[(scores_df['tissue'] == 'Leaf') & (scores_df['treatment'].isin(['Control', 'FeLim']))]
            scores_df = scores_df[['rxn_ID', 'day', 'treatment', 'bind', 'norm_value',  \
                        'rxn_score_I_dist', 'subsystems', 'rxn_dist_quantile', 'features']]

            scores_df['spc'] = spc
            scores_df['res'] = msr
            scores_df['spc_res'] = spc+'_'+msr
            scores_df['res_spc'] = msr+'_'+spc
            # print("----------------------------- "+spc+'_'+msr)
            if all_df.empty:
                 all_df = scores_df
            else: 
                all_df = pa.concat([all_df, scores_df], ignore_index=True)
            # print(all_df.shape)

    all_rxn_scores(RPHelper, all_df)












