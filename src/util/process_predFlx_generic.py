import sys
from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from src.util.parameters import *
import src.util.process_RES_generic as prg
import src.reaction_scores.reactionScoresHelper as rsh

import pandas as pa
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import plotly.io as pio
import seaborn as sns
import plotly.graph_objects as go

A, B = 1, -1

time_points = ["02d", "04d", "07d", "14d", "21d"]
# time_points = ["21d"]
treatments = ["Control", "FeLim", "FeEX", "ZnLim", "ZnEx"]
tissue = "Leaf"
# spc = "Sorghum"
spc = "Sorghum"
trmt_colm = 'treatment'
value_col = 'value'
sv = 15
plot = False

omic_mm_results = f"/Users/selalaoui/Projects/AMN/omic_amn_mm/Result/sv{sv}"
figures_path =  f"/Users/selalaoui/Projects/AMN/omic_amn_mm/Result/pred_fluxes_figures/sv{sv}/"
             # /Users/selalaoui/Projects/AMN/omic_amn_mm/Result

def get_rxn_ID(row):
    if any(y in row['rxn_ID'] for y in ['_f', '_r', '_i', '_o']):
        # print(row['rxn_ID'].rsplit("_", 1)[0])
        id_only = row['rxn_ID'].rsplit("_", 1)[0]
    else:
        id_only = row['rxn_ID']

    return id_only

def plot_rxn_scores_subsystem(rxn_scores_df, group_cols=['treatment', 'day', 'tissue'], trmt_col='treatment', control_name='Control', score_perct=90, project='QPSI', verbose=False):
    subsystem_check = False
    rxn_scores_df.drop(columns=['value'], inplace=True)
    rxn_scores_df.rename(columns={'norm_value':'value'}, inplace=True)
    rxn_scores_df['tissue'] = 'Leaf'

    group_cols.remove(trmt_col)
    cols_indices = [i for i in range(0, len(group_cols))]
    print(group_cols)
    # cols_indices.remove(group_cols.index(trmt_col))
    print(cols_indices)

    dist = np.abs(rxn_scores_df['rxn_score_I_dist'])\
        .describe([score_perct/100])[str(score_perct)+'%']
    dist_dict = dict()
    # group_cols = [group_cols[idx] for idx in cols_indices]
    # print(group_cols)
    groups = rxn_scores_df.groupby(group_cols)
    print(groups)

    for name, group in groups:
        name = (name[0], name[0]) if name[-1] == '' else (name[0], name[-1])
        dist_dict[name] = np.abs(group['rxn_score_I_dist'])\
            .describe([score_perct/100])[str(score_perct)+'%']
        dist_dict[name] = dist
    print(dist_dict)

    # print(abc)
    if verbose: print(f"Distance at {score_perct} score_perct: {dist}")
    print(f"Distance at {score_perct} score_perct: {dist}")

    if subsystem_check:
        rxn_scores_df['subsystems'] = rxn_scores_df.apply(lambda row: apply_literal_eval(row), axis=1)
        rxn_scores_df = rxn_scores_df.explode('subsystems')
        sys_class = get_subsysClass(list(rxn_scores_df['subsystems'].unique()))
        ## remove
        subsystems = []
        if not subsystems:
            subsystems = ["Core_tetrapyrrole_biosynthesis_in_plants",
            "Chlorophyll_Biosynthesis_in_plants_and_prokaryotes"]#,
        print(subsystems)
        print(rxn_scores_df.shape)

        # rxn_scores_df = rxn_scores_df[rxn_scores_df['subsystems'].isin(subsystems)]
        rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants_and_prokaryotes','')
        rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_in_plants','')
        rxn_scores_df['subsystems'] = rxn_scores_df['subsystems'].str.replace('_',' ')
        # # if you need to remove some subsystems
        # remove_sub = [...]
        # rxn_scores_df = rxn_scores_df[~rxn_scores_df['subsystems'].isin(remove_sub)]
        rxn_scores_df['class'] = rxn_scores_df['subsystems'].map(sys_class)

    rxn_scores_df[(np.abs(rxn_scores_df['rxn_score_I_dist']) >=dist)].to_csv(f"{spc}_{score_perct}Perc_reactions.csv", index=False)
    print(rxn_scores_df[(np.abs(rxn_scores_df['rxn_score_I_dist']) >=dist)]['rxn_ID'].unique())
    # print(abc)

    # Change from long to wide
    col = [trmt_col]
    val = ['value', 'rxn_score_I_dist', 'rxn_dist_quantile']#,
    ind = set(rxn_scores_df.columns) - set().union(col, val)
    distance_df = rxn_scores_df.pivot(index=ind, columns=col, values=val)
    distance_df.columns = distance_df.columns.get_level_values(0) + '_' +  distance_df.columns.get_level_values(1)
    distance_df = distance_df.reset_index()

    if verbose: print(distance_df.columns)


    if project == 'hAlpha':
        wt, ht = 1400, 1000
    elif project == 'QPSI':
        wt, ht = 1400, 400
    else:
        wt, ht = 1650, 600

    names = {'value_'+trmt : trmt for trmt in treatments}
    names[f'value_{control_name}'] = control_name
    distance_df.rename(columns=names, inplace=True)
    if 'time_stamp' in  distance_df:
        distance_df.sort_values(by=['time_stamp'], inplace=True)
    else:
        distance_df.sort_values(by=['day'], inplace=True)

    for trmt in treatments:
        # Normalized scores scatter plots
        mx = distance_df["rxn_score_I_dist_"+trmt].max()
        mx = max(np.abs(distance_df["rxn_score_I_dist_"+trmt].min()), mx)
        # mx = 0.25
        title = f"{spc} {trmt} -- Normalized predicted fluxes ({score_perct} percentile)."

        category_orders = {} if project == "secMeta" \
                            else {"day": time_points} # tissues}

        fig = px.scatter(
                distance_df,
                # x=trmt,
                # y=control_name,
                x=control_name,
                y=trmt,
                color="rxn_score_I_dist_"+trmt,
                title=title,
                # symbol= 'subsystems',
                color_continuous_scale='icefire',
                range_color=[-1*mx, mx],
                facet_col="day",
                # facet_row=trmt_col,
                hover_data=["rxn_ID"],
                labels={"rxn_score_I_dist_"+trmt: "Flux Dist"},
                category_orders=category_orders
                , height=ht, width=wt
                # , facet_row_spacing=0.08
                , facet_col_spacing=0.03
            )
        reference_line_x_range = np.array([0, 1])


        col_values = time_points
        row_values = "Leaf" # treatments
        for row_idx, row_figs in enumerate(fig._grid_ref):
            for col_idx, col_fig in enumerate(row_figs):
                # Add the greater variability subsystems reactions
                col_val = col_values[col_idx] if isinstance(col_values, list) else col_values
                row_val = row_values[(row_idx+1)*-1] if isinstance(row_values, list) else row_values

                trmt_df = distance_df[(np.abs(distance_df['rxn_score_I_dist_'+trmt]) >= \
                                        dist_dict[(col_val, row_val)])
                                      & (distance_df[group_cols[cols_indices[0]]] == col_val)
                                      & (distance_df[group_cols[cols_indices[-1]]] == row_val)]

                fig = fig.add_trace(go.Scatter(x=trmt_df[control_name],
                                y=trmt_df[trmt],
                                # x=sys_trmt_1[trmt],
                                # y=sys_trmt_1[control_name],
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

        # fig.update_xaxes(range=[-0.01, 1.05], showticklabels=True)
        # fig.update_yaxes(range=[-0.01, 1.05])#, showticklabels=True)

        fig.update_xaxes(range=[-0.01, 0.6], showticklabels=True)
        fig.update_yaxes(range=[-0.01, 0.6])#, showticklabels=True)

        # fig.update_xaxes(range=[-0.01,0.5], showticklabels=True)
        # fig.update_yaxes(range=[-0.01,0.5], showticklabels=True)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=12
            ), showlegend=False
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
         # \
         #                                        if len(a.text.split("="))>1)\
         #                                        else a.text.split("=")[0]))
        # fig.for_each_trace(lambda t: t.update(name=t.name.split("=")[1] \
        #                                     if len(t.name.split("="))>1 \
        #                                     else t.name.split("=")[0]))
        # fig.show()
        # fig.update_xaxes(title_text=treatment, showticklabels=True, zerolinewidth=1, zerolinecolor='#000000')
        # fig.update_yaxes(title_text=control_name, zerolinewidth=1, zerolinecolor='#000000')

        # print(len(fig.data))
        ln_rows = len(row_values) if isinstance(row_values, list) else 1
        num_cells = len(col_values) + ln_rows
        fig.data = [fig.data[i] for i in reversed(range(len(fig.data)))]

        fig.show()
        # plot_path = omic_mm_results+f"/scatter_plots_{trmt}.png"
        plot_path = figures_path+f"{spc}_scatter_plots_{trmt}.png"
        pio.write_image(fig, plot_path, scale=2, width=wt, height=ht)

    compute_var_std(rxn_scores_df, trmt_col)
    # rxn_scores_df.rename(columns={'value':'norm_value'})

def compute_var_std(input_df, trmt_col, cprmt='', var_column='rxn_score_I_dist', verbose=False):
    print(" --------- Processing RIPTiDe results * Fluxes scatter std line plots --------- ")

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
    sort_dict = {'tissue': ['Leaf', 'Root'], 'treatment': ['FeEX', 'FeLim', 'ZnEx', 'ZnLim'], 'day': ['02d', '04d', '07d', '14d', '21d']}
    group_cols = ['treatment', 'day'] #['tissue', 'growth_stage', 'treatment']


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

    # up_trans = up_std.T
    # up_trans.columns = up_trans.columns.get_level_values(0)
    up_trans = up_std[sort_dict[group_cols[-1]]]
    up_trans = up_trans.T
    print(up_trans.head())
    # print(abc)


    if verbose: print(f"downregulated std")
    down_std = df[df[var_column] <= 0][group_cols+[var_column]].\
        groupby(group_cols).var().unstack(1)
    down_std.columns = down_std.columns.get_level_values(1)
    # down_trans = down_std.T
    # down_trans.columns = down_trans.columns.get_level_values(0)
    down_trans = down_std[sort_dict[group_cols[-1]]]
    down_trans = down_trans.T
    print("-----------> up")
    print(down_trans.head())
    # print(abc)

    upmax = up_trans.select_dtypes(include=[np.number]).max().max()
    downmax = down_trans.select_dtypes(include=[np.number]).max().max()
    ymax = max(upmax, downmax)
    if verbose: print(f"---------------  Maximum STD {ymax}")

    print(up_trans.head())
    print(group_cols)
    for trmt in sort_dict[group_cols[0]]:
        col = f'{trmt}'
        if verbose: print(f"Generating figure for {col}")
        width, height = 4.5, 2
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
        ax.set_ylim([0, ymax+0.000005])
        # ax.legend(frameon=True, loc='best', facecolor='green', framealpha=1, fancybox=True)
        plt.legend(frameon=True, loc='best', facecolor='white', framealpha=.8, fancybox=True)
        fig.suptitle(f"{col.replace('_', ' ')} {cprmt}")

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.savefig(figures_path+f'{spc}_std_var_{col}.png', bbox_inches='tight', dpi=400)
        # plt.show()
    #
    # for tissue in sort_dict[group_cols[-1]]:
    #     for gs in sort_dict[group_cols[0]]:
    #         col = f'{gs}_{tissue}'
    #         if verbose: print(f"Generating figure for {col}")
    #         width, height = 4.5, 2
    #         fig, ax = plt.subplots(figsize=(width, height))
    #         ax.plot(up_trans[col], label='Up-regulation', color='chocolate')
    #         ax.plot(down_trans[col], label='Down-regulation', color='steelblue')
    #         # ax.legend(fancybox=True)
    #         # leg=plt.legend(loc='best', numpoints=1, fancybox=True)
    #         ax.grid()
    #         ax.ticklabel_format(style='sci', axis='y', useOffset=True)
    #         # plt.ticklabel_format(style='sci', useOffset=True)
    #         plt.xlabel('Treatment')
    #         plt.ylabel(r'$\sigma$')
    #         ax.set_ylim([0, ymax])
    #         # ax.legend(frameon=True, loc='best', facecolor='green', framealpha=1, fancybox=True)
    #         plt.legend(frameon=True, loc='best', facecolor='white', framealpha=.8, fancybox=True)
    #         fig.suptitle(f"{col.replace('_', ' ')} {cprmt}")
    #
    #         ax.spines['top'].set_visible(False)
    #         ax.spines['right'].set_visible(False)
    #         ax.spines['bottom'].set_visible(False)
    #         ax.spines['left'].set_visible(False)
    #         plt.savefig(f'std_figures/{col}.png', bbox_inches='tight', dpi=400)
    #         # plt.show()


fluxes_df = pa.DataFrame()
for tp in time_points:
    more = "_noADP_noP"
    tp_df = pa.read_csv(f"{omic_mm_results}/{spc}_{tissue}_{tp}_complexFix{more}_V_rxn_fba_Vbf_RES.tsv", sep="\t")
    tp_df = tp_df[["rxn_ID", "treatment", "Pred", "subsystems"]]
    tp_df.rename(columns={"Pred": value_col}, inplace=True)
    tp_df['day'] = tp

    group_cols = [trmt_colm, 'day']
    tp_df = rsh.compute_rxn_variability(tp_df, group_cols, treatments, treatments[0], trmt_colm, value_col)

    fluxes_df = pa.concat([fluxes_df, tp_df], ignore_index=True)

print(fluxes_df.head())
# print(abc)

 
if plot:
    score_perct=98
    plot_rxn_scores_subsystem(fluxes_df, group_cols=['treatment', 'day', 'tissue'], trmt_col='treatment', control_name='Control', score_perct=score_perct, project='QPSI', verbose=False)


#####  reaction scores
res_path = project_root+"integration_results/reaction_scores_binding_Jul2/"\
+f"plastidial_model/{spc}_objective_abundance_Control.tsv"
res_df = pa.read_csv(res_path, sep="\t")
res_df = res_df[(res_df['tissue'] == tissue) & (res_df['time_stamp'].isin(time_points))]
res_df = res_df[["treatment", "time_stamp", "value", "subsystems", "rxn_ID", "norm_value", "rxn_score_I_dist", "rxn_dist_quantile"]]
res_df = res_df.rename(columns = {"value": "score", "time_stamp": "day", "norm_value": "score_norm", "rxn_dist_quantile": "score_dist_quantile"})
print(res_df.head())

#####  reaction fluxes
print(fluxes_df.columns)
# print(abc)
fluxes_df.rename(columns = {"value": "pred_flux", "norm_value": "pred_flux_norm", "rxn_score_I_dist": "flux_I_dist",  "rxn_dist_quantile": "flux_dist_quantile"}, inplace=True)
fluxes_df['rxn_ID_only'] = fluxes_df.apply(lambda row: get_rxn_ID(row), axis=1)
print(fluxes_df.head())
print(fluxes_df.columns)

##### merge dataframes
fluxes_df = res_df.merge(fluxes_df[["treatment", "day", "rxn_ID_only", "rxn_ID", "pred_flux",
                                    "pred_flux_norm", "flux_I_dist",  "flux_dist_quantile"]],
                        how="outer",
                        left_on=["treatment", "day", "rxn_ID"],
                        right_on=["treatment", "day", "rxn_ID_only"],
                        suffixes=("_x", None))
fluxes_df.drop(columns=['rxn_ID_only', "rxn_ID_x"], inplace=True)
print(fluxes_df.head())
fluxes_df[["rxn_ID", "treatment", "day", "score", "score_norm", "rxn_score_I_dist", "score_dist_quantile", "pred_flux", "pred_flux_norm", "flux_I_dist", "flux_dist_quantile", "subsystems"]].to_csv(f"{spc}_leaf_res_flux.tsv", sep='\t', index=False)
