import cobra as co
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.flux_analysis import variability as va
from cobra.io import read_sbml_model

import pandas as pa

## check if a reaction is a transport reaction
def is_transport(row, co_model):
    if 'EX_' in row.rxn_ID:
        return True

    try:
        rxn = co_model.reactions.get_by_id(row.rxn_ID)
        if len(rxn.compartments) > 1:
            return True
        else:
            return False

    except KeyError:
        print("reaction {} not found".format(row.rxn_ID))
        return False

## print model stats 
def stats(model, tissue):
    ex_num = 0
    d_trans_num = 0
    c_trans_num = 0
    plast_num = 0
    cyto_num = 0
    for rxn in model.reactions:
        if 'EX_' in rxn.id:
            ex_num += 1
        elif len(rxn.compartments) > 1:
            if 'd0' in rxn.id: d_trans_num+=1
            else: c_trans_num +=1
            # trans_num += 1
        else:
            if 'd0' in rxn.id: plast_num += 1
            else: cyto_num += 1


    print("-*"*10 + "-  "+tissue)
    print(f"There are:")
    print(f"  {len(model.metabolites)} metabolies")
    print(f"  {len(model.genes)} genes")
    print(f'  {len(model.medium)} media')
    print(f"  {len(model.reactions)} reactions")
    print("   Reaction Stats:")
    print(f"    {ex_num} export reactions")
    print(f"    {d_trans_num+c_trans_num} transporters {d_trans_num} plastidial and {c_trans_num} cyto")
    print(f"    {plast_num} plastidial reactions")
    print(f"    {cyto_num} cyto reactions")

## Compute FVA and FBA  results for the metabolic model
## returns: pandas dataframe with reaction ID, reaction flux, FVA maximum and minimum fluxes,
##             and if the reaction is a transport reaction.
def run_FVA(xmlModel_path, tissue, spc='', metabolites=False):
    fva_df = pa.DataFrame()

    # Read COBRApy model
    co_model = read_sbml_model(xmlModel_path)
    # print model stats
    stats(co_model, tissue)

    # Get reaction flux from FBA simulation
    solution = co_model.optimize()
    solution = solution.fluxes.to_frame()
    solution[solution.abs() < 10**-6] = 0
    solution = solution.rename_axis('rxn_ID').reset_index()

    # find the list of active metabolites
    if metabolites:
        all_mets = set()
        rxns = solution[solution['fluxes'] > 0]['rxn_ID'].unique()
        for rxn_id in rxns:
            rxn = co_model.reactions.get_by_id(rxn_id)
            all_mets.update(["{}, {}".format(met.id.split('_')[0], met.name) for met in rxn.metabolites])
        print(f"{spc}_{tissue}_active_metabolites -- > {len(all_mets)}")
        with open(f'{spc}_{tissue}_active_metabolites.txt', 'w') as f:
            for met in all_mets:
                f.write(f"{met}\n")

    # Run FVA
    fva_df = fva(co_model, co_model.reactions, processes=1)
    if fva_df.empty:
        return fva_df
    # Set very small fluxes to 0
    fva_df[fva_df.abs() < 10**-6] = 0
    # Add isTransport column
    fva_df = fva_df.rename_axis('rxn_ID').reset_index()
    fva_df['isTrans'] = fva_df.apply(lambda row: is_transport(row, co_model), axis=1)
    # Add flux to the dataframe
    fva_df =  pa.merge(fva_df, solution, how="left", on='rxn_ID')
    fva_df['tissue'] = tissue
    fva_df = fva_df[['rxn_ID', 'fluxes', 'maximum', 'minimum', 'isTrans', 'tissue']]
    return fva_df
