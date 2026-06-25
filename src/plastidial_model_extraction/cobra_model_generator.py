import sys
from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)

import src.util.modelComponents as mc
from src.util.modelComponents import bcolors

import cobra as co
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.flux_analysis import variability as va
import pandas as pa
import json

class CobraHelper:
    def __init__(self, plast_model, co_model=None):
        if isinstance(plast_model, mc.Model):
            self.plast_model = plast_model
        else:
            raise ValueError(bcolors.FAIL+" Wrong model type provided."+bcolors.ENDC)

        if co_model == None: co_model = co.Model()
        self.co_model = co_model

        self.solution = None
        self.blocked = set()
        self.co_blocked = set()

    ## --- Print model details (reactions, genes, etc.) and summarized model properties
    #      (number of reactions, genes, etc.)
    def print_details(self):
        print("Reactions")
        print("---------")
        for x in self.co_model.reactions:
            print("%s : %s" % (x.id, x.reaction))

        print("")
        print("Metabolites")
        print("-----------")
        for x in self.co_model.metabolites:
            print('%9s : %s' % (x.id, x.formula))

        print("")
        print("Genes")
        print("-----")
        for x in self.co_model.genes:
            associated_ids = (i.id for i in x.reactions)
            print("%s is associated with reactions: %s" %
                  (x.id, "{" + ", ".join(associated_ids) + "}"))

    def print_properties(self):
        print(bcolors.PROG)
        print(f'\t {len(self.co_model.reactions)} reactions')
        print(f'\t {len(self.co_model.metabolites)} metabolites')
        print(f'\t {len(self.co_model.genes)} genes')
        print(f'\t {len(self.co_model.medium)} media')
        print(f'\t {len(self.co_model.compartments)} compartments')
        print(bcolors.ENDC)
        # print(self.co_model.objective)
    ## ----------------------------------------------------------------------------------------

    # convert JSON str reversibility to COBRApy upper and lower bounds.
    @staticmethod
    def get_reversibility(direction):
        config = co.Configuration()

        if direction == '<':
            return config.lower_bound, 0
        elif direction == '>':
            return 0, config.upper_bound
        else: # '='
            return config.lower_bound, config.upper_bound

    # Generate SBML GPR rules based on the JSON "modelReactionProteins" -- Recursive Method
    def generate_gpr(self, reactionProteins, key_level):
        operand_dict={
            "modelReactionProteins": " or ",
            "modelReactionProteinSubunits": " and ",
            "feature_refs": " or "
        }
        next_key_dict = {
            "modelReactionProteins": "modelReactionProteinSubunits",
            "modelReactionProteinSubunits": "feature_refs"
        }

        if isinstance(reactionProteins, str):
            return reactionProteins.split('/')[-1]

        if isinstance(reactionProteins, mc.Feature):
            return reactionProteins.gene_id

        if isinstance(reactionProteins, list):
            if len(reactionProteins) == 0:
                return ''

            if key_level in next_key_dict:
                next_key = next_key_dict[key_level]
                element_list = [self.generate_gpr(getattr(pe, next_key), next_key) for pe in reactionProteins]
            else:
                element_list = [self.generate_gpr(pe, "") for pe in reactionProteins]
            # clean the list: remove None elements
            element_list = [e for e in element_list if (isinstance(e, (str, list)) and len(e)>0)]
            rule = operand_dict[key_level].join(element_list)

            if len(element_list) > 1 : rule = '(' + rule + ')'
            return rule

        return ''

    # Add all reactions, genes, metabolites, etc. to the new COBRApy model
    def create_cobra_model(self, media_from_json, media_cprmt, cprmts):
        co_model = co.Model()
        rxn_list = dict()
        metabolites_dict = dict()

        # ADD metabolites
        print(bcolors.SUBRESULT+"    Adding metabollites")
        for cpd in self.plast_model.modelcompounds_dict.values():
            metabolites_dict[cpd.id] = co.Metabolite(
                                cpd.id,
                                formula=cpd.formula,
                                name=cpd.name,
                                charge=cpd.charge,
                                compartment=cpd.compartment)

        # ADD media
        print(bcolors.SUBRESULT+"    Adding media exchange reactions")
        for cpd in media_from_json["mediacompounds"]:
            # 1. read media and create metabolite
            if cpd["id"]+media_cprmt not in metabolites_dict:
                new_cpd = None
                for cpt in cprmts:
                    try:
                        new_cpd = self.plast_model.modelcompounds_dict[cpd["id"]+cpt].deepcopy()
                        new_cpd.update_cpd_compartment(media_cprmt)
                        break
                    except KeyError:
                        continue

                if new_cpd != None:
                    print("adding new metabolites {}".format(new_cpd.id))
                    metabolites_dict[cpd["id"]+media_cprmt] = co.Metabolite(
                            new_cpd.id,
                            formula=new_cpd.formula,
                            name=new_cpd.name,
                            charge=new_cpd.charge,
                            compartment=new_cpd.compartment)
                else:
                    print("adding new metabolites {}".format(cpd["id"]+media_cprmt))
                    metabolites_dict[cpd["id"]+media_cprmt] = co.Metabolite(
                            cpd["id"]+media_cprmt,
                            formula=cpd["name"],
                            name=cpd["name"],
                            compartment=media_cprmt.replace('_', ''))
            # 2. create the reaction
            rxn_id = 'EX_'+cpd["id"]+media_cprmt
            reaction = co.Reaction(id = rxn_id)
            reaction.name = "Exchange for "+cpd["name"]
            # reaction.lower_bound , reaction.upper_bound = get_reversibility('=')
            reaction.lower_bound = int(cpd["minFlux"]) if cpd["minFlux"] > -1000 else co.Configuration().lower_bound
            reaction.upper_bound = int(cpd["maxFlux"]) if cpd["minFlux"] < 1000 else co.Configuration().upper_bound
            # 3. add the metabolite to the reaction
            reaction.add_metabolites({metabolites_dict[cpd["id"]+media_cprmt] : -1})

            rxn_list[rxn_id] = reaction


        # ADD reactions
        print(bcolors.SUBRESULT+"    Adding reactions")

        no_ands, no_ors, no_singles, no_none = (0, 0, 0, 0)
        for rxn in self.plast_model.modelreactions_dict.values():

            reaction = co.Reaction(id = rxn.id)
            reaction.name = rxn.name
            reaction.lower_bound , reaction.upper_bound = self.get_reversibility(rxn.direction)
            # Append all metabolites
            rxn_metabolites_dict = dict()
            for cpd in rxn.modelReactionReagents:
                cpd_id = cpd.id
                rxn_metabolites_dict[metabolites_dict[cpd_id]] = cpd.coefficient
            reaction.add_metabolites(rxn_metabolites_dict)
            # Append all genes
            if rxn.modelReactionProteins != None:
                # # If COBRApy allows complex rules
                rule = self.generate_gpr(getattr(rxn, 'modelReactionProteins'), "modelReactionProteins")

                if 'and' in rule.lower():  no_ands += 1
                elif 'or' in rule.lower(): no_ors += 1
                elif rule == '':          no_none += 1
                else:                     no_singles += 1

                if rule != "":
                    reaction.gene_reaction_rule = rule

            rxn_list[rxn.id] = reaction

        sum = no_ands+no_ors+no_singles+no_none


        # Add objective reaction
        print(bcolors.SUBRESULT+"    Adding biomass reaction")
        obj_id = 'bio1_biomass'
        if obj_id not in self.plast_model.biomasses_dict:
            biomass = self.plast_model.biomasses_dict['bio1']
        else:
            biomass = self.plast_model.biomasses_dict[obj_id]
        obj = co.Reaction(id = obj_id)
        obj.name = biomass.name
        obj.lower_bound , obj.upper_bound = self.get_reversibility('>')

        # Append all metabolites
        obj_metabolites_dict = dict()
        for cpd in biomass.biomasscompounds:
            cpd_id = cpd.id
            obj_metabolites_dict[metabolites_dict[cpd_id]] = cpd.coefficient
        obj.add_metabolites(obj_metabolites_dict)
        rxn_list[obj_id] = obj

        co_model.add_reactions(list(rxn_list.values()))
        co_model.objective = obj_id

        print(bcolors.SUBRESULT+"    Finalizing model"+bcolors.ENDC)
        self.co_model = co_model
        self.print_properties()
        bcolors.ENDC

        return co_model

    # Run COBRApy FBA optimization to check if the model is functional
    def run_optimization(self):
        self.solution = self.co_model.optimize()
        print(bcolors.RESULT+"  {}".format(self.solution)+bcolors.ENDC)

    # Write reaction fluxes to file:
    #   rxn_id rxn_definition flux max_flux min_flux
    def write_flux_details(self, output_folder:str=''):
        fluxes_df = self.co_model.optimize().fluxes
        fva_df = fva(self.co_model, self.co_model.reactions, processes=1)
        # Set reactions fluxes < 10^-6 to 0
        # fluxes_df[fluxes_df.abs() < 10**-6] = 0
        # print(fluxes_df.head())
        # fva_df[fva_df.abs() < 10**-6] = 0

        print(bcolors.PROMPT+" Writing fluxes details to "+output_folder+"model_fluxes.tsv"+bcolors.ENDC)
        with open(output_folder+'model_fluxes.tsv','w') as fh:
            fh.write("rxn_id\tdefinition\tflux\tmax_flux\tmin_flux\n")
            for i in range(len(fluxes_df.index)):
                rxn_cpt = fluxes_df.index[i]
                rxn = None
                cpd = None
                if('cpd' in rxn_cpt):
                    cpd = rxn_cpt.split('_', 1)[1]
                else:
                    rxn = rxn_cpt

                flux = "{:.2f}".format(fluxes_df.iloc[i])

                rxn_fva = fva_df.loc[fva_df.index == rxn_cpt]
                min = 'N/A'
                max = 'N/A'
                if not rxn_fva.empty:
                    min = "{:.2f}".format(rxn_fva.iloc[0]['minimum'])
                    max = "{:.2f}".format(rxn_fva.iloc[0]['maximum'])

                if(rxn is not None):
                    model_reaction = self.co_model.reactions.get_by_id(rxn).build_reaction_string(use_metabolite_names=True)
                    fh.write(rxn_cpt+"\t"+str(model_reaction)+"\t"+str(flux)+"\t"+str(max)+"\t"+str(min)+"\n")
                elif(cpd is not None):
                    fh.write(rxn_cpt+"\t"+self.co_model.metabolites.get_by_id(cpd).name+"\t"+str(flux)+"\t"+str(max)+"\t"+str(min)+"\n")
                else:
                    fh.write(rxn_cpt+"\t"+str(flux)+"\t"+str(max)+"\t"+str(min)+"\n")

    # Write list of biomass metabolites, their names and coefficients to file
    def write_biomass_details(self, output_folder=''):
        from BiochemPy import Compounds
        cpd_helper = Compounds()

        print(bcolors.PROMPT+" Writing biomass details to "+output_folder+"model_biomass.tsv"+bcolors.ENDC)
        bio=self.co_model.reactions.get_by_id("bio1_biomass")
        with open(output_folder+'model_biomass.tsv','w') as fh:
            fh.write("Metabolite\tname\tcoefficient\n")
            for metabolite, coeff in bio.metabolites.items():
                fh.write('%s \t %s \t %s \n' % (metabolite, metabolite.name, coeff))

    # Read KBase FBA results from file, run COBRApy "find_blocked_reactions" to find blocked
    # reactions using both methods.
    def FBA_analysis(self, KBase_FBA, verbose=False):
        # Reactions blocked on COBRApy FBA
        print(bcolors.PROMPT+" Running COBRApy FVA"+bcolors.ENDC)
        COBRA_blocked = list()
        COBRA_blocked = va.find_blocked_reactions(self.co_model, self.co_model.reactions, open_exchanges=True, processes=1)
        print(bcolors.RESULT+" {} reactions blocked in COBRApy simulation.".format(len(COBRA_blocked)))
        if verbose: print(COBRA_blocked)

        # Reactions blocked on KBase FBA
        kbase_fba = pa.read_csv(KBase_FBA, sep='\t')
        KBase_blocked = kbase_fba[(kbase_fba['min_flux']==0) & (kbase_fba['max_flux']==0)]
        ## REMOVE
        KBase_blocked = kbase_fba[(kbase_fba['min_flux']>1000) & (kbase_fba['max_flux']==0)]
        print(" {} reactions blocked in KBase FBA.".format(KBase_blocked.shape[0]))

        # reactions blocked in KBase and COBRApy
        blocked_overlap = KBase_blocked[KBase_blocked['id'].isin(COBRA_blocked)]
        blocked_overlap_list = blocked_overlap['id'].tolist()
        print(" {} common blocked reactions.".format(blocked_overlap.shape[0]))

        # blocked in COBRApy and not in KBase
        cobra_extra = set(COBRA_blocked) - set(blocked_overlap_list)
        print(" {} extra reactions blocked in COBRApy.".format(len(cobra_extra)))

        common_trans = set()
        for rxn_id in blocked_overlap_list:
            rxn = self.co_model.reactions.get_by_id(rxn_id)
            if len(rxn.compartments) > 1:
                common_trans.add(rxn_id)
        if verbose: 
            print(bcolors.SUBRESULT+"   Common blocked transporters: {}".format(common_trans))

        extra_trans = set()
        for rxn_id in cobra_extra:
            rxn = self.co_model.reactions.get_by_id(rxn_id)
            if len(rxn.compartments) > 1:
                extra_trans.add(rxn_id)

        print("   Blocked transporters in COBRApy only {}".format(extra_trans)+bcolors.ENDC)

        self.blocked.update(common_trans)
        self.co_blocked.update(extra_trans)
