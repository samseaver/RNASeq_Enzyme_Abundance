import sys

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from src.util.modelComponents import *

from cobra_model_generator import CobraHelper

from cobra.io import write_sbml_model
import pandas as pa
import json
import os

import BiochemPy

output_folder = project_root+'integration_results/plastidial_models_details/'
param_file = project_root+"src/plastidial_model_extraction/parameters.json"

printModelDetails = False
writeModelXML = False
writeModelAnalysis = False

# Helper class that extracts the plastidial model from the full JSON model
# it read parameters and reaction information from the parameters.json file,
# which must be modified depending on model/compartment specific requirements.
class ModelBuilder:
    def __init__(self, json_file:str=None, xml_file:str=None):
        if(json_file is None and xml_file is None):
            raise ValueError(bcolors.FAIL+"Model Builder requires the path to a ModelSEED JSON or a COBRA SBML file."+bcolors.ENDC)
        
        self.model_json = None
        self.model_xml = None
        self.output_model_file = None
        if(json_file):
            self.model_json = json_file
            self.output_model_file = "plastidial-"+os.path.splitext(json_file)[0]
            p = Path(json_file)
            self.output_model_file = str(p.parent / f"plastidial-{p.stem}")

        if(xml_file):
            self.model_xml = xml_file
            p = Path(xml_file)
            self.output_model_file = str(p.parent / f"plastidial-{p.stem}")

        self.missing_rxns_dict = dict()
        self.main_cprmt = None
        self.extra_cprmts = list()

        self.media_cprmt = None
        self.cprmts = None
        self.missing_trans = None
        self.remove_trans = None
        self.reactionsFromMSD = None
        self.exclude_reactions = []

        self.folder = None        
        
        self.media_json = None
        self.load_parameters_from_file(param_file)

        self.full_model = None
        # Load the full model from file
        print(bcolors.PROMPT+" Reading full model JSON:"+bcolors.ENDC)
        print("\t"+bcolors.PROG+self.model_json+bcolors.ENDC)
        with open(self.model_json, 'r') as f:
            data=f.read()
            json_model = json.loads(data)
            self.full_model = Model().fromJSON(json_model)
        if self.full_model == None:
            raise ValueError(bcolors.FAIL+"  Full model couldn't load."+bcolors.ENDC)

        # Init the plastidial model to be built
        self.plast_model = None
        self.media_from_json = None

        self.plast_reactions = list()
        self.transport_plast_rxns_list = list()

        self.plast_products_set = set()
        self.plast_reactants_set =set()
        self.extra_products_set = set()
        self.extra_reactants_set =set()

        self.trans_for = set()
        self.media_for = dict()

    ## Load project-specific parameters from file
    #   and create output folders
    def load_parameters_from_file(self, param_file):
        param_dict = None
        with open(param_file, 'r') as f:
            data=f.read()
            param_dict = json.loads(data)
        if param_dict == None:
            raise ValueError(bcolors.FAIL+"  Couldn't load parameters."+bcolors.ENDC)

        # Compartments
        self.main_cprmt = param_dict["compartments"]["main_cprmt"]
        self.extra_cprmts = param_dict["compartments"]["extra_cprmts"]
        self.media_cprmt = param_dict["compartments"]["media_cprmt"]
        self.cprmts = param_dict["compartments"]["cprmts"]

        # Reactions
        self.missing_trans = param_dict["reactions"]["missing_trans"]
        self.remove_trans = set(param_dict["reactions"]["remove_trans"])
        self.reactionsFromMSD = set(param_dict["reactions"]["reactionsFromMSD"])
        self.missing_rxns_dict = dict()
        rxn_set = set()
        for subsys in param_dict["reactions"]["main_cprmt_rxns"]:
            rxn_set.update(param_dict["reactions"][subsys])
        self.missing_rxns_dict[self.main_cprmt] = rxn_set

        for cprt, subsys_list in param_dict["reactions"]["extra_cprmt_rxns"].items():
            rxn_set = set()
            for subsys in subsys_list:
                rxn_set.update(param_dict["reactions"][subsys])
            # {rxn for subsys_list in subsys
            #         for rxn in param_dict["reactions"][subsys]}
            self.missing_rxns_dict[cprt] = rxn_set

        # Biomass
        # self.biomass_file = project_root+param_dict["files_paths"]["biomass_file"]

    def load_media_file(self, media_file:str=None):
        with open(media_file, 'r') as f:
            data=f.read()
            self.media_from_json = json.loads(data)
        if self.media_from_json == None:
            raise ValueError(bcolors.FAIL+"  Couldn't load media from file."+bcolors.ENDC)

    # extract all the plastidial, exchange, sink, demand, transport and biomass reactions
    def filter_reactions(self, printModelDetails=False):
        for rxn_id, rxn in self.full_model.modelreactions_dict.items():
            if rxn_id in self.exclude_reactions:
                continue

            rxnType = rxn.isTransporter()

            if rxnType == 'OtherTransport' and rxn.compartment in self.main_cprmt:
                self.transport_plast_rxns_list.append(rxn)
                # Add non-plastidial reagents to the trans_for set
                for id in rxn.reagents:
                    self.trans_for.add(id)
                    
                    # Store internal transporters in media_for as well!
                    # (You'd need to strip the compartment to match the base metabolite ID)
                    base_id = id.split('_')[0] 
                    media_id = base_id + self.media_cprmt
                    self.media_for.setdefault(media_id, set()).add(rxn)
                continue
            
            if rxnType == 'MediumTransport' and (rxn.compartment in [self.main_cprmt.replace('_', '')]+[cprt.replace('_', '') for cprt in self.extra_cprmts]):
                # add media reactants
                for id in rxn.reagents:
                    if id.endswith(self.media_cprmt): 
                        # Initialize a set if missing, then add the reaction
                        self.media_for.setdefault(id, set()).add(rxn)
                continue

            if rxn.compartment in self.main_cprmt:
                self.plast_reactions.append(rxn)
                self.plast_reactants_set.update(rxn.reactants)
                self.plast_products_set.update(rxn.products)
                continue

        if printModelDetails:
            print(bcolors.PROG+"\tPlastidial transport reactions: {}".format(len(self.transport_plast_rxns_list)))
            # print("\tMedia transport reactions: {}".format(len(self.transport_media_rxns_list)))
            print("\tPlastidial reactions: {}.".format(len(self.plast_reactions))+bcolors.ENDC)

    @staticmethod
    def compound_from_MSD(cpd_id):
        from BiochemPy import Compounds
        cpd_helper = Compounds()
        cpd_dict = cpd_helper.loadCompounds()

        try:
            cpd = cpd_dict[cpd_id.split('_')[0]]
            return Compound(id=cpd_id, name=cpd["abbreviation"],
                    formula=cpd["formula"], charge=cpd["charge"])
        except KeyError:
            print(bcolors.WARNING+"\t Couldn't find compound {} in MSD -- Skipping".format(cpd_id)+bcolors.ENDC)

        return None

    # adds reaction to the new model
    #  -> Add reagents from full model to new plastidial model
    #       if reagent is not found, looks it up in the ModelSEEDDB
    #  -> Adds compartment to the plastid model
    def add_reaction_toModel(self, rxn):
        if rxn.id in self.plast_model.modelreactions_dict:
            return

        new_rxn = rxn.deepcopy()
        self.plast_model.add_new_reaction(new_rxn)

        # Add reaction metabolites to the plastidial model
        for m_id in new_rxn.reagents:
            if m_id not in self.plast_model.modelcompounds_dict:
                try:
                    new_cpd = self.full_model.modelcompounds_dict[m_id].deepcopy()
                except KeyError:
                    new_cpd = self.compound_from_MSD(m_id)

                self.plast_model.add_new_metabolite(new_cpd)

        # Add reaction compartment to the plastidial model
        for cprmt in new_rxn.compartments:
            if cprmt not in self.plast_model.modelcompartments_dict:
                try:
                    new_cpt = self.full_model.modelcompartments_dict[cprmt].deepcopy()
                except KeyError:
                    new_cpt = Compartment(cpmt_letterId=cprmt[:-1],
                    cpmt_index=cprmt[-1], label="Unknown")
                self.plast_model.add_new_compartment(new_cpt)

    # Add full model plastidial reaction to new model
    #   -> List generated using "filter_reactions"
    #   -> Uses "add_reaction_toModel" to add the reactions
    def add_plastidial_reactions(self):
        for rxn in self.plast_reactions:
            self.add_reaction_toModel(rxn)
        print(bcolors.PROG+"\t {} plastidial reaction(s) added to model ({})".format(len(self.plast_reactions), self.main_cprmt.replace('_', ''))+bcolors.ENDC)

    # Reactions from other compartments that are needed for the model to work
    #  are added here.
    def add_missing_reactions(self, cprmt, extra=True, verbose=False):
        still_missing = self.missing_rxns_dict[cprmt]
        num_added = 0
        other_cprmts = self.cprmts.copy()
        if verbose:
            print(extra)
            print("Processing : ", cprmt, other_cprmts)

        if (cprmt in other_cprmts):
            other_cprmts.remove(cprmt)

        if 'y' in cprmt:
            other_cprmts = [cprmt]

        if verbose: print(other_cprmts)

        for o_cprmt in other_cprmts:
            added = set()
            for rxn_id in still_missing:
                if extra and verbose:
                    print("Processing ", rxn_id, " ", o_cprmt)

                # if reaction is already in the model, add it directly
                if rxn_id+cprmt in self.full_model.modelreactions_dict:
                    if extra and verbose:
                        print("Adding ", rxn_id, " ", o_cprmt)

                    added.add(rxn_id)
                    num_added += 1
                    self.add_reaction_toModel(self.full_model.modelreactions_dict[rxn_id+cprmt].deepcopy())
                    continue

                # Try finding the reaction in o_cprmt
                try:
                    new_rxn = self.full_model.modelreactions_dict[rxn_id+o_cprmt].deepcopy()
                    isTrans = new_rxn.isTransporter()
                    new_rxn.update_reaction_compartment(cprmt, isTrans)
                    # add any missing metabolites to the model
                    for m_id in new_rxn.reagents:
                        self.full_model.add_metabolite_toCprmt(m_id, o_cprmt, cprmt)

                    added.add(rxn_id)
                    num_added += 1
                    self.add_reaction_toModel(new_rxn)
                except KeyError:
                    continue

            still_missing = still_missing - added

        if still_missing:
            print(bcolors.PROG+'\tReactions not in the full model: {}'.format(still_missing)+bcolors.ENDC)

        print(bcolors.PROG+"\t {} missing reaction(s) added to model ({})".format(num_added, cprmt.replace('_', ''))+bcolors.ENDC)

    # Reactions not included in the full model, but that are needed in the plastidial model
    #  can be included from the ModelSEEDDB
    def add_reactions_form_ModelSEEDDB(self, cpmt):

        # fetch reactions from MDB
        print(bcolors.PROG+"\tFetching reactions from MSD")
        reactions_helper = BiochemPy.Reactions()
        reactions_dict = reactions_helper.loadReactions()
        num_added = 0
        for rxn_id in self.reactionsFromMSD:
            reaction = None
            try:
                reaction =  reactions_dict[rxn_id]
            except KeyError:
                print(bcolors.PROG+"\t "+reaction["id"]+" not found in MSD")

            if reaction != None:
                print(bcolors.PROG+"\t Creating a new entry for "+rxn_id+cpmt)
                num_added += 1
                # ADD reagents
                cpd_list = reactions_helper.parseEquation(reaction['code'])
                r_list = list()
                for cpd in cpd_list:
                    cpd_id = cpd["compound"]+cpmt
                    r_list.append(Reagent(cpd["coefficient"],cpd_id))
                    # add metabolite to model if not in
                    if not (cpd_id in self.full_model.modelcompounds_dict):
                        print(bcolors.PROG+"\t    Adding missing metabolite: "+cpd_id)
                        m = Compound(id=cpd_id, name=cpd["name"],
                                formula=cpd["formula"], charge=cpd["charge"])
                        self.full_model.add_new_metabolite(m)
                # create new reaction
                new_reaction_id = reaction["id"]+cpmt
                new_reaction = Reaction(id=new_reaction_id, direction=reaction["reversibility"],
                modelReactionReagents=r_list, name=reaction["abbreviation"],
                probability=0, protons=0)

                self.add_reaction_toModel(new_reaction)


        print(bcolors.PROG+"\t {} reaction(s) added to model ({}) from MSD".format(num_added, cpmt.replace('_', ''))+bcolors.ENDC)

    # Add new transport reactions that are needed for the model
    # either for new exchange reactions (see add_media_exchange_reactions)
    # or non-plastidial transporters (see add_transport_reactions)
    def new_transport(self, cpd_id, cpd_name, cprmt1, cprmt2):
        rxn_id = 'rxn08217' if cpd_id == 'cpd00137' else "rxn_"+cpd_id

        r_media = Reagent(-1, cpd_id+cprmt1)
        r_extra = Reagent(1, cpd_id+cprmt2)

        new_trans = Reaction(id=rxn_id+cprmt2, direction='=',
        modelReactionReagents=[r_media, r_extra],
        modelReactionProteins=None, name="Transport from "+cpd_name)
        # self.plast_model.add_new_reaction(new_trans)

        self.trans_for.add(cpd_id+cprmt1)
        self.trans_for.add(cpd_id+cprmt2)
        return new_trans

    # read the media compounds from the media file
    #   add exchange reactions for all plastid media compounds
    def add_media_exchange_reactions(self, extra_cprmt):
        print(bcolors.SUBRESULT+"    Adding (media) exchange reactions"+bcolors.ENDC)
        num_added = 0
        for cpd in self.media_from_json["mediacompounds"]:
            
            # ==========================================================
            # PRE-SEED METABOLITES: Ensure the compound exists in ALL 3 compartments
            # ==========================================================
            for cprt in [self.media_cprmt, extra_cprmt, self.main_cprmt]:
                cpd_id_cprt = cpd["id"] + cprt
                if cpd_id_cprt not in self.plast_model.modelcompounds_dict:
                    new_cpd = None
                    for search_cprt in [self.media_cprmt, self.main_cprmt] + self.extra_cprmts:
                        search_id = cpd["id"] + search_cprt
                        if search_id in self.full_model.modelcompounds_dict:
                            new_cpd = self.full_model.modelcompounds_dict[search_id].deepcopy()
                            break
                    if new_cpd == None:
                        new_cpd = self.compound_from_MSD(cpd["id"])
                    if new_cpd != None:
                        new_cpd.update_cpd_compartment(cprt)
                        self.plast_model.add_new_metabolite(new_cpd)

            # ==========================================================
            # EVALUATE NATIVE EXCHANGE REACTIONS (Now handling sets!)
            # ==========================================================
            has_media_exchange = False
            has_plastid_transport = False
            
            if cpd["id"]+self.media_cprmt in self.media_for:
                m_rxns = self.media_for[cpd["id"]+self.media_cprmt] # <--- Retrieves the set!
                
                for m_rxn in m_rxns:
                    m_rxn_cprts = [c.replace('_', '') for c in m_rxn.compartments]

                    # Media exchange reaction goes to extra_cprmt (e.g. c0)
                    if extra_cprmt.replace('_', '') in m_rxn_cprts:
                        self.add_reaction_toModel(m_rxn)
                        num_added += 1
                        has_media_exchange = True
                        
                        # Check if we need to build the downstream Cytosol <=> Plastid leg
                        if cpd["id"]+self.main_cprmt not in self.trans_for:
                            new_rxn = self.new_transport(cpd["id"], cpd["name"], extra_cprmt, self.main_cprmt)
                            self.transport_plast_rxns_list.append(new_rxn)
                            self.trans_for.add(cpd["id"]+self.main_cprmt)
                            has_plastid_transport = True

                    # Media exchange reaction goes directly to main_cprmt (e.g. d0)
                    elif self.main_cprmt.replace('_', '') in m_rxn_cprts:
                        self.add_reaction_toModel(m_rxn)
                        num_added += 1
                        has_media_exchange = True
                        has_plastid_transport = True # Skips the middleman entirely!
                        
                        if cpd["id"]+extra_cprmt not in self.trans_for:
                            new_rxn = self.new_transport(cpd["id"], cpd["name"], extra_cprmt, self.main_cprmt)
                            self.transport_plast_rxns_list.append(new_rxn)
                            self.trans_for.add(cpd["id"]+extra_cprmt)

                    # The exchange reaction goes to a totally different compartment
                    else:
                        for cpt in [extra_cprmt, self.main_cprmt]:
                            new_rxn = self.new_transport(cpd["id"], cpd["name"], self.media_cprmt, cpt)
                            self.transport_plast_rxns_list.append(new_rxn)
                            for m_id in new_rxn.reagents:
                                self.trans_for.add(m_id)
                        has_media_exchange = True
                        has_plastid_transport = True

            # ==========================================================
            # ADD COMPLETELY SYNTHETIC EXCHANGES IF NOTHING WAS FOUND
            # ==========================================================
            if not has_media_exchange:
                print(bcolors.PROG+"\t   Adding new exchange reaction for "+cpd["id"]+bcolors.ENDC)
                m_rxn_e = self.new_transport(cpd["id"], cpd["name"], self.media_cprmt, extra_cprmt)
                self.transport_plast_rxns_list.append(m_rxn_e)

            # If the inner leg is still missing after evaluating all native reactions
            if not has_plastid_transport and cpd["id"]+self.main_cprmt not in self.trans_for:
                print(bcolors.PROG+"\t   Adding new transport reaction for "+cpd["id"]+bcolors.ENDC)
                m_rxn_c = self.new_transport(cpd["id"], cpd["name"], extra_cprmt, self.main_cprmt)
                self.transport_plast_rxns_list.append(m_rxn_c)
                self.trans_for.add(cpd["id"]+self.main_cprmt)

        print(bcolors.PROG+"\t {} media transport reaction(s) added to model".format(num_added)+bcolors.ENDC)

    # Add alist of missing transporters provided in the parameters file: reactions->missing_trans
    def add_transport_reactions(self, verbose=False):
        num_added = 0
        num_skip = 0
        # Add missing reactions
        for cprmt, cpd_ids in self.missing_trans.items():
            for cpd_id in cpd_ids:
                rxn = self.new_transport(cpd_id, cpd_id, self.main_cprmt, cprmt)
                self.transport_plast_rxns_list.append(rxn)

        # Add model reaction
        for rxn in self.transport_plast_rxns_list:
            if rxn.idOnly in self.remove_trans:
                print(bcolors.PROG+'\tSkipping transport reaction: '+rxn.id)
                num_skip += 1
                continue

            if (rxn.id == "rxn10832_d0") and verbose:
                print("\tadding rxn10832_d0 to model")

            # check if all metabolites are in model
            reagents = rxn.reagents
            if all(cpd in self.plast_model.modelcompounds_dict for cpd in reagents):
                self.add_reaction_toModel(rxn)
                num_added += 1
            else:
                for id in reagents:
                    self.trans_for.discard(id)
                # print(reagents)
                # print(bcolors.WARNING+"\t Skipping transporter {}: ".format(rxn.id)+
                #     "not all reagents are in model."+bcolors.ENDC)
                num_skip += 1

        print(bcolors.PROG+"\t {} transporters added to model".format(num_added)+bcolors.ENDC)
        print(bcolors.PROG+"\t {} transporters skipped bacause of missing compounds".format(num_skip)+bcolors.ENDC)

    # Add the biomass reaction:
    #   - read plastid biomass compounds from file defined in the parameters
    #               file (files_paths->biomass_file)
    #   - read all other information from the full model biomass reactions
    def add_biomass_bio1(self, biomass_file, biomass_id:str='bio1'):
        bio_full_model = self.full_model.biomasses_dict[biomass_id]

        print(bcolors.PROMPT+" Loading plastidial biomass ..."+bcolors.ENDC)
        plast_biomass_df = pa.read_csv(biomass_file)

        plast_biomass_df['id_d'] = plast_biomass_df['id'].astype(str)+'_'+plast_biomass_df['compartment']+'0'
        plast_biomass_dict = dict(zip(plast_biomass_df.id_d, plast_biomass_df.stoichiometry))

        cpd_list = list()
        for cpd_id, coeff in plast_biomass_dict.items():
            if cpd_id in self.plast_model.modelcompounds_dict:
                cpd_list.append(Reagent(coeff=float(coeff)*-1, cpd_id = cpd_id))
            else:
                print(bcolors.PROG+"\tBiomass cpd not in model: {}".format(cpd_id)+bcolors.ENDC)

        new_bio = Biomass(bio_full_model.id, "Plant Plastid", cpd_list,
        bio_full_model.cellwall, bio_full_model.cofactor, bio_full_model.dna,
        bio_full_model.energy, bio_full_model.lipid, bio_full_model.other,
        bio_full_model.protein, bio_full_model.rna)

        self.plast_model.add_new_biomass(new_bio)

    # Write the model to JSON file
    def write_json_file(self, fileName="new_model.json"):
        json_output = json.dumps(self.plast_model.as_dict(), indent=4)
        output_json_file = fileName
        with open(output_json_file, 'w') as outfile:
            outfile.write(json_output)

    # Call reactions above to add different types of reactions, build the new model
    #   and write it to file, if write2json is True.
    def construct_new_model(self, write2json=False):
        self.plast_model =  Model(name=None, id=None, modelReactions=dict(),
                        biomasses=dict(), modelCompartments=dict(),
                        modelCompounds=dict(), source=None, source_id=None,
                        gapfillings=None, gapgens=None, template_ref=None,
                        genome_ref=None, type=None)
        # Get all reactions (compounds and compartments)
        print(bcolors.SUBRESULT+"    Adding full model {} reactions".format(self.main_cprmt)+bcolors.ENDC)
        self.add_plastidial_reactions()
        # for key, rxn in self.plast_model.modelreactions_dict.items():
        #     print(rxn)
        print(bcolors.SUBRESULT+"    Adding missing {} reactions".format(self.main_cprmt)+bcolors.ENDC)
        self.add_missing_reactions(self.main_cprmt, extra=False)
        for cpmt in self.extra_cprmts:
            print(bcolors.SUBRESULT+"    Adding extra missing reactions to {}".format(cpmt)+bcolors.ENDC)
            self.add_missing_reactions(cpmt, extra=True)
        print(bcolors.SUBRESULT+"    Adding {} reactions from MSD".format(self.extra_cprmts[0])+bcolors.ENDC)
        self.add_reactions_form_ModelSEEDDB(self.extra_cprmts[0])
        print(bcolors.SUBRESULT+"    Adding transport reactions"+bcolors.ENDC)
        self.add_transport_reactions()

        print(bcolors.SUBRESULT+"    Finalizing model ..."+bcolors.ENDC)
        # Updating model attributes
        self.plast_model.set_model_attributes(self.full_model.name,
                self.full_model.id, self.full_model.source, self.full_model.source_id, self.full_model.gapfillings, self.full_model.gapgens, self.full_model.template_ref, self.full_model.genome_ref, self.full_model.type)

        if write2json: self.write_json_file()

    # Write model statistics (biomass atoms fractional sums) to file after
    #  building the model
    def write_biomass_atoms_details(self, output_folder=''):
        print(bcolors.PROMPT+" Writing biomass atoms fractional sums "+output_folder+"biomass_atom_fraction.tsv"+bcolors.ENDC)
        # loop over full model biomass compounds
        print(bcolors.PROG+"    Extracting full biomass compounds and computing fractions"+bcolors.ENDC)
        full_atom_dict = self.full_model.compute_atoms_franctions(self.biomass_id)

        print(bcolors.PROG+"    Extracting plastidial biomass compounds and computing fractions"+bcolors.ENDC)
        plastid_atom_dict = self.plast_model.compute_atoms_franctions(self.biomass_id)

        print(bcolors.PROG+"    Writing fractions to file"+bcolors.ENDC)
        with open(output_folder+'biomass_atom_fraction.tsv','w') as fh:
            fh.write("atom\tbiomass\tplastid_biomass\t%\n")
            all_atoms = set().union(full_atom_dict.keys(), plastid_atom_dict.keys())

            for atom in all_atoms:
                full_sum = 0 if atom not in full_atom_dict else full_atom_dict[atom]
                plastid_sum = 0 if atom not in plastid_atom_dict else plastid_atom_dict[atom]
                percent = 'N/A' if full_sum == 0 else (plastid_sum* 100) / full_sum
                fh.write('%s \t %.2f \t %.2f \t %.2f \n' % (atom, full_sum, plastid_sum, percent))

# Helper class that runs a pipline to build
#    1. the JSON plastidial model
#    2. the SBML model using COBRApy
#    3. clean both models from inactive reactions (FVA min and max fluxes == 0)
#    4. write both models to file
class ModelGenerator:
    def __init__(self, modelbuilderobject):
        self.mb = modelbuilderobject
        self.cobrahelperobject = None

    # Utilize ModelBuilder class to read full model, extract reactions of interest
    #  and build the new model
    def generate_model(self):
        # Get the list of reactions
        print(bcolors.PROMPT+" Filtering model reactions "+bcolors.ENDC)
        self.mb.filter_reactions()
        # Construct the new JSON model
        print(bcolors.PROMPT+" Constructing plastidial model "+bcolors.ENDC)
        self.mb.construct_new_model()
        print(bcolors.PROG+"\t  reactions {}".format(len(self.mb.plast_model.modelreactions_dict)))
        print(bcolors.PROG+"\t  Compounds {}".format(len(self.mb.plast_model.modelcompounds_dict)))
        print(bcolors.PROG+"\t  Compartment {}".format(len(self.mb.plast_model.modelcompartments_dict))+bcolors.ENDC)

    # Build SBML model using the JSON model build in the previous step using
    #  use class "CobraHelper" from "cobra_model_generator.py"
    def cobra_model(self):
        # Create COBRApy model
        print(bcolors.PROMPT+" Converting JSON model to COBRA model"+bcolors.ENDC)
        # Create a COBRA helper object
        self.cobrahelperobject = CobraHelper(self.mb.plast_model)

        self.cobrahelperobject.create_cobra_model(self.mb.media_from_json,
        self.mb.media_cprmt,[self.mb.main_cprmt]+self.mb.extra_cprmts)

        print(bcolors.PROMPT+" Running optimization"+bcolors.ENDC)
        # self.cobrahelperobject.run_optimization()
        # self.cobrahelperobject.FBA_analysis(self.mb.KBase_FBA)

        # print model properties and details
        if printModelDetails :
            self.cobrahelperobject.print_properties()
            self.cobrahelperobject.print_details()

    # Use FVA results to remove blocked reactions from the model (i.e. min and max fluxes == 0)
    #  then write both model to JSON and xml files.
    def clean_write_model(self, clean_up:bool=False):
        if clean_up and self.cobrahelperobject is not None:
            print(bcolors.PROMPT+" Cleaning up the models "+bcolors.ENDC)
            num_removed = 0
            print(bcolors.PROG+"\tremoving the following blocked reactions: ", self.cobrahelperobject.co_blocked, bcolors.ENDC)
            for rxn_id in self.cobrahelperobject.blocked:
                # if 'y' in rxn_id: continue
                try:
                    self.cobrahelperobject.co_model.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)
                    del self.mb.plast_model.modelreactions_dict[rxn_id]
                    num_removed += 1
                except KeyError:
                    continue

            for rxn_id in self.cobrahelperobject.co_blocked:
                # if 'y' in rxn_id: continue
                self.cobrahelperobject.co_model.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)
                num_removed += 1
                del self.mb.plast_model.modelreactions_dict[rxn_id]
                num_removed += 1

            print(bcolors.PROG+"\t Removed {} blocked reaction(s) across the two models.".format(num_removed)+bcolors.ENDC)

            print(bcolors.PROMPT+" Running optimization after cleaning model"+bcolors.ENDC)
            self.cobrahelperobject.run_optimization()

        # write model to file
        print(bcolors.PROMPT+" Writing models to file:"+bcolors.ENDC)
        json_model = self.mb.plast_model.as_dict()
        with open(self.mb.output_model_file+".json", "w") as outfile:
            print("\t"+self.mb.output_model_file+".json")
            json.dump(json_model, outfile, indent=4)

        # Save it as xml
        if writeModelXML and self.cobrahelperobject is not None:
            write_sbml_model(self.cobrahelperobject.co_model,
            self.mb.output_model_file+".xml")

        if(self.cobrahelperobject is not None):
            self.cobrahelperobject.print_properties()

        # Print model details:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        if writeModelAnalysis:
            self.mb.write_biomass_atoms_details(output_folder)
            if(self.cobrahelperobject is not None):
                self.cobrahelperobject.write_biomass_details(output_folder)
                self.cobrahelperobject.write_flux_details(output_folder)

    # Run the pipline
    def run_model_generator(self, toJSON = True, toXML = False, clean_up = False):
        self.generate_model()

        if toXML:
            self.cobra_model()

        if toJSON:
            self.clean_write_model(clean_up)

if __name__ == '__main__':
    mBuilder = ModelBuilder()
    mGen = ModelGenerator(mBuilder)
    mGen.run_model_generator(toFile = True)
