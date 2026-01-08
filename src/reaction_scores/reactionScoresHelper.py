import warnings
warnings.simplefilter(action='ignore', category=Warning)

import pandas as pa
import sys
import numpy as np

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from src.util.bcolors import bcolors

# Compute the protein subunit score using one of three methods:
#   1. relab -- relative abundance: computes the fraction of each paralog and multiply
#                                   it by its weight: weight_frac. All paralogs weight_frac
#                                   values are summed up to return the subunit score.
#   2. Sum: subunit score is the sum of the abundances of its paralogs
#   3. Max: subunit score is the max of the abundances of its paralogs
def compute_subunit_score(df, id_col, value_col, std_col, group_cols, method='relab', rxn_id = ''):

    rxn_id_only = rxn_id.split('_')[0]

    if method == 'relab':
        # df['weight'] = df['me'].div(df.groupby(['date', 'rank'])['me'].transform('sum'))
        df['relab_frac'] = df[value_col].div(df.groupby(group_cols)[value_col].transform('sum'))

        df['weight_frac'] = df['relab_frac']*df['weight']

        
        sums_df =  df.groupby(group_cols).agg({id_col:','.join, \
            'weight_frac':'sum', value_col: 'sum'\
            }).reset_index()


    elif method == 'sum':
    
        sums_df =  df.groupby(group_cols).agg({value_col: 'sum' \
                        , id_col:','.join \
                        # , id_col: lambda x: ','.join(str(x))
                        # , id_col: lambda x: list(x) \lambda x: ','.join(x))
                        ,std_col: lambda std_col : np.sqrt((std_col*std_col).sum()),\
                    }).reset_index()
    else:
        sums_df = df.loc[df.groupby(group_cols)[value_col].idxmax()]

    return sums_df


# Compute the protein complex's score using one of two methods:
#   1. relab -- relative abundance: multiply each subunit weight by the minimum relative abundance
#                                   across all subunits return the sums representing the protein
#                                   weight
#   2. Min: subunit score is the min of its sub-units
def compute_mrp_score(df, value_col, group_cols, method='relab', rxn_id = ''):
    if method == 'relab':
        # multiply each subunit weight by the minimum relative abundance across all subunits
        # return the sums representing the protein weight
        sums_df =  df.groupby(group_cols).apply(lambda df,weight,relab: \
            sum(df[weight] * df[relab].min()), 'weight_frac', value_col)
        sums_df = sums_df.reset_index(name='weightedSumWeight')

        rxn_id_only = rxn_id.split('_')[0]

    else:
        sums_df = df.loc[df.groupby(group_cols)[value_col].idxmin()]

    return sums_df


# Compute the protein complex's score using one of two methods:
#   1. relab -- relative abundance: sum the scores (weightedSumWeight) of all reacrion protein
#                                   complex
#   2. Sum: sum the scores of all reacrion protein complex
#   3. Max: max of all reacrion protein complex scores
# Add features, binding information and subsystem data to the result

def genera_list(rows):
    return [y for x in rows
                for y in x.split(',')
            ]


def compute_rxn_score_value(df, id_col, rxn, value_col, std_col, group_cols, subsystems_set, ftr_list, binding, method='relab', rxn_id = ''):
    new_df = pa.DataFrame()

    if method == 'relab':
        new_df = df.groupby(group_cols).agg({'weightedSumWeight':'sum'}).reset_index()
        
    elif method == 'sum':
        new_df = df.groupby(group_cols).agg({value_col:'sum'\
                , id_col: ','.join \
                # , id_col: lambda x : genera_list(x) \
                ,std_col: lambda std_col : np.sqrt((std_col*std_col).sum()),\
                }).reset_index()
    else:
        new_df = df.loc[df.groupby(group_cols)[value_col].idxmax()]

    # new_df.rename(columns = {value_col:'value'}, inplace=True)
    new_df = new_df.assign(subsystems = [subsystems_set for i in new_df.index])

    if not binding:
        binding = "None"
    new_df = new_df.assign(bind = [binding for i in new_df.index])
    new_df['rxn_ID'] = rxn
    new_df['features'] = str(ftr_list)
    # new_df["bind"] = binding
    rxn_id_only = rxn_id.split('_')[0]

    return new_df


# Itertes over all model reactions and computes reaction score using:
#   - compute_subunit_score
#   - compute_mrp_score
#   - compute_rxn_score_value
def compute_model_score(relab_data, metModel, id_col, value_col, std_col, group_cols, spc_name, method='relab', verbose=False):
    # -- input
    # relab_data.columns: gene_id, treatment, tissue, time_stamp, value (TPM, TMM, or any other)
    # model_rxn_dict: rxn_id -> Reaction object (modelComponents)
    # -- output
    # dfcolumns: rxn_ID, treatment, subsystems, 2d, 4d, ...
    nr = 0
    rxn_scores = pa.DataFrame()
    n_rxn = len(metModel.modelreactions_dict)

    # Implementing EXTREAM algorithm here
    # Create a set of unique (Gene, Base_Reaction) tuples
    # We strip the last 2 chars if they are directional tags (_f or _r) to get the biological ID
    # Using a set automatically handles the deduplication (counting the pair only once)
    unique_gene_rxn_pairs = set(
        (g, r.id[:-2] if r.id.endswith(('_f', '_r')) else r.id)
        for r in metModel.modelreactions_dict.values()
        for g in r.genes
    )

    # Count reaction associations per gene
    gene_rxns_count = pa.Series([pair[0] for pair in unique_gene_rxn_pairs]).value_counts()

    # Store original values and Normalize
    relab_data['original_value'] = relab_data[value_col]
    relab_data[value_col] /= relab_data[id_col].map(gene_rxns_count).fillna(1)

    for rxn_id, rxn in metModel.modelreactions_dict.items():
        nr+=1
        rxn_subsys = rxn.subsystems

        if not rxn.genes:
            # print(bcolors.PROG+" -- No features, skipping {} ({}/{})".format(rxn_id, nr, n_rxn)+bcolors.ENDC)
            continue

        # print(bcolors.PROG+"Computing scores for {} ({}/{})".format(rxn_id, nr, n_rxn)+bcolors.ENDC)

        mdlrxn_prot_list = rxn.gpr

        ftr_list = list()
        mrp_scores = pa.DataFrame() # <-- return max
        # for plotting: 
        prot = 0 
        rxn_id_only = rxn_id.split('_')[0]

        for mdlrxn_prot in mdlrxn_prot_list:
            prot += 1
            mrps_scores = pa.DataFrame() # <-- return min
            sub = 0
            for mdlrxn_subunit_ftrs in mdlrxn_prot:
                sub += 1

                mdlrxn_subunit_ftrs = [ftr for ftr in mdlrxn_subunit_ftrs if not any(pref in ftr.upper() for pref in ['ATC', 'ATM'])]
                
                if not mdlrxn_subunit_ftrs:
                    continue
    
                feature_value_df = relab_data[(relab_data[id_col].isin(mdlrxn_subunit_ftrs))]
                feature_value_df[id_col] = feature_value_df[id_col].astype(str)
                ftr_list.extend(mdlrxn_subunit_ftrs)
                mrps_score = compute_subunit_score(feature_value_df, id_col, value_col, std_col, group_cols, method, rxn_id = f"{rxn_id}_{spc_name}_{prot}_{sub}")
                mrps_scores= pa.concat([mrps_scores, mrps_score], ignore_index=True)

            if not mrps_scores.empty:
                mrp_score  = compute_mrp_score(mrps_scores, value_col, group_cols, method, rxn_id = f"{rxn_id}_{spc_name}_{prot}")
                mrp_scores= pa.concat([mrp_scores, mrp_score], ignore_index=True)

        if not mrp_scores.empty:
            rxn_score = compute_rxn_score_value(mrp_scores, id_col, rxn_id, value_col, std_col, \
                                        group_cols, sorted(rxn_subsys), ftr_list, rxn.binds, method, rxn_id = rxn_id)
            
            rxn_scores= pa.concat([rxn_scores, rxn_score], ignore_index=True)

    if verbose: print(rxn_scores.head(5))
    return rxn_scores


# For every treatment, compute the distance between the identity line and the point:
#       x <- treatment score
#       y <- control score
def compute_rxn_variability(scores_df, group_cols, treatments, control_id, trmt_colm, value_col, percentile=90, verbose=False):
    if verbose: print(" --------- Computing rxn score varibility --------- norm_"+value_col)
    quantiles = [i*5/100 for i in range(0, 21)]
    labels = [int(i*100) for i in quantiles[:-1]]
    # Identity line
    A, B = 1, -1
    result = pa.DataFrame()

    # Normalize:
        # create groups list
    gcs = group_cols.copy()
    gcs.remove(trmt_colm)
    scores_df[value_col] = scores_df[value_col].astype('float')
        # Group the scores DF
    groups = scores_df.groupby(gcs)
        # computes group-wise mean/std,
        # then auto broadcasts to size of group chunk
    min = groups[value_col].transform("min")
    max = groups[value_col].transform("max")
    scores_df['norm_'+value_col] = (scores_df[value_col] - min) / (max - min)

    if verbose: print(scores_df.head(5))
    # Compute distance
    for trmt in treatments:
        trmtScores_df = scores_df[scores_df[trmt_colm].isin([trmt, control_id])].copy()
        trmtScores_df = trmtScores_df[['norm_'+value_col, 'rxn_ID']+group_cols]

        #    pivot the DF
        col = [trmt_colm]
        val = ['norm_'+value_col]
        ind = set(trmtScores_df.columns) - set().union(col, val)
        trmtScores_df = trmtScores_df.pivot(index=ind, columns=col, values=val)

        #    reset reset columns names
        trmtScores_df.columns = trmtScores_df.columns.get_level_values(0) + '_' +  trmtScores_df.columns.get_level_values(1)
        trmtScores_df = trmtScores_df.reset_index()

        #    Compute distance
        trmtScores_df['rxn_score_I_dist'] = np.abs(A*trmtScores_df['norm_'+value_col+'_'+trmt] + B*trmtScores_df['norm_'+value_col+'_'+control_id]) / np.sqrt(np.power(A,2)+np.power(B,2))
        #    Set the direction of the change
        trmtScores_df.loc[(trmtScores_df['norm_'+value_col+'_'+trmt] < trmtScores_df['norm_'+value_col+'_'+control_id]), "rxn_score_I_dist"] = trmtScores_df["rxn_score_I_dist"]*-1

        #    Add treatment column
        trmtScores_df[trmt_colm] = trmt

        #    Append to full DF
        result = pa.concat([result, trmtScores_df[['rxn_score_I_dist', trmt_colm]+list(ind)]], ignore_index=True)
        if verbose: print(f"Appending {result.columns}")

    scores_df = pa.merge(scores_df, result, how='left', on=['rxn_ID']+group_cols)

    # ### #### Compute quantiles
    scores_df['rxn_dist_quantile'] = 0
    scores_df['rxn_score_I_dist_abs'] = np.abs(scores_df['rxn_score_I_dist'])
    scores_df['rxn_dist_quantile'] = scores_df['rxn_score_I_dist_abs']\
        .rank(method='first')\
        .transform(lambda x: pa.qcut(np.abs(x), quantiles, labels=labels))

    if verbose: print(scores_df.head(5))
    scores_df.drop(["rxn_score_I_dist_abs"], axis=1, inplace=True)
    if verbose: print(scores_df.columns)

    return scores_df
