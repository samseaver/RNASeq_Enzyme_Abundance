import sys
import os
import csv
import json
from urllib.request import urlopen
import numpy as np
import pandas as pa

import cobra as co
from cobra.io import read_sbml_model, write_sbml_model
from time import time
# from cobra.summary.metabolite_summary import MetaboliteSummary as ms


from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
import src.util.modelComponents as mc
from src.util.bcolors import bcolors

class bcolors:
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

class ModelDataProcessingHelper:
    def __init__(self, species='ptrich', ortho_version='apr10'):
        self.species = species

        # self.models_folder = os.path.join(project_root, "data", "metabolic_models", "plastidial_models")
        self.models_folder = os.path.join(project_root, "data", "metabolic_models", "fullmodels_media")
        self.version = ortho_version 
        self.orthologs_folder = os.path.join(project_root, "data", "orthologs", f"orthologs_{self.version}")

        self.model_name_json = "none"
        self.atha_model_json = 'none'
        self.orthologs_file  = "none"

        self.co_model = None
        self.spc_json_model = dict()
        self.atha_json_model = dict()
        self.allFeatures = set()
        self.ortho_mapping_dict = dict()

        self.ortho_dict = dict()

        fileNames = os.listdir(self.models_folder)
        for fileName in fileNames:
            if all((y in fileName.lower()) for y in [self.species.lower(), 'json', 'ComplexFix'.lower()]) :#and ('fds' in fileName): #, 'sucrose'
                self.model_name_json = fileName
                print(self.model_name_json)
            if all((y in fileName.lower()) for y in ['atha', 'json', 'ComplexFix'.lower()]): #and ('fds' in fileName):
                self.atha_model_json = fileName
                print(self.atha_model_json)

            if (self.model_name_json != 'none') and (self.atha_model_json != 'none'):
                break

        if (self.model_name_json == 'none') or (self.atha_model_json == 'none'):
            print(bcolors.WARNING+f"Couldn't find metabolic models, exiting!"+bcolors.ENDC)
            sys.exit(1)

    def readModelJSON(self):
        print(bcolors.PROG+f"Reading model ..."+bcolors.ENDC)
        # other species
        self.spc_json_model = dict()
        model_json_file_path = os.path.join(self.models_folder, self.model_name_json)
        with open(model_json_file_path, 'r') as f:
            self.spc_json_model = mc.Model.fromJSON(json.load(f))

        # Arabidopsis model
        model_json_file_path = os.path.join(self.models_folder, self.atha_model_json)
        with open(model_json_file_path, 'r') as f:
            self.atha_json_model = mc.Model.fromJSON(json.load(f))

        if not (self.spc_json_model and self.spc_json_model):
            print(bcolors.WARNING+"Unable to read model, exiting!"+bcolors.ENDC)
            sys.exit(1)

    def readOrthologs(self, mode=1, version='v'):
        # Try to fing the orthologs files in the folder 'self.orthologs_folder' defined in the init method
        self.orthologs_file = 'none'
        fileNames = os.listdir(self.orthologs_folder)
        print(fileNames)
        for fileName in fileNames:
            if (self.species.lower() in fileName.lower()) \
                and ('homolog' in fileName.lower()) \
                and (version in fileName.lower()):
                self.orthologs_file = fileName
                break
        
        if self.orthologs_file == 'none':
            print(bcolors.WARNING+f"No orthologs file found for {self.species.upper()}, exiting!"+bcolors.ENDC)
            sys.exit(1)

        print(bcolors.PROG+f"Generating orthologs dict from {self.orthologs_file} ..."+bcolors.ENDC)
        self.orthologs_file = os.path.join(self.orthologs_folder, self.orthologs_file)
        sep = '\t'

        with open(self.orthologs_file, mode='r') as inFile:
            reader = csv.reader(inFile, delimiter=sep)
            for row in reader:
                if (mode==1) and (float(row[5]) <= 0.5):
                    continue

                if (mode==2) and (row[6] != 'O'):
                    continue

                if (mode==3) and (row[6] not in ['O', 'PAH', 'PAO', 'PAM', 'PA']):
                    continue 


                if 'at' in row[3].lower():
                    AT_gene = row[3].split('.')[0]
                    spc_gene = row[4].rsplit('.', 2)[0]
                else:
                    AT_gene = row[4].split('.')[0]
                    spc_gene = row[3].rsplit('.', 2)[0]

                self.ortho_dict[spc_gene] = (row[0], row[6])

                # # skipping organellar genes
                if any(pref in AT_gene.upper() for pref in ['ATC', 'ATM']):
                    continue
                try:
                    self.ortho_mapping_dict[AT_gene].add(spc_gene)
                except KeyError:
                    self.ortho_mapping_dict[AT_gene] = {spc_gene}


        if not self.ortho_mapping_dict:
            print(bcolors.WARNING+"No ortholog mapping found!"+bcolors.ENDC)
            sys.exit(1)

        print(f"Retrieved {len(self.ortho_mapping_dict)} orthologs ...")
        print(list(self.ortho_mapping_dict.items())[:5])

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

        if isinstance(reactionProteins, list):
            if len(reactionProteins) == 0:
                return ''

            if key_level in next_key_dict:
                next_key = next_key_dict[key_level]
                element_list = [self.generate_gpr(pe[next_key], next_key) for pe in reactionProteins]
            else:
                element_list = [self.generate_gpr(pe, "") for pe in reactionProteins]
            # clean the list: remove None elements
            element_list = [e for e in element_list if (isinstance(e, (str, list)) and len(e)>0)]
            rule = operand_dict[key_level].join(element_list)

            if len(element_list) > 1 : rule = '(' + rule + ')'
            return rule

        return ''

    def replaceOrthologsInModel(self, vebose=True):
        rxn_summary_str = ""
        newGenes, notfound = set(), set()
        allModelGenes, allRemoved, allNew = set(), set(), set()
        no_pred_in_new = set()
        no_pred_in_old = set()
        no_pred_at_all = set()
        rxns_atha_pred, spc_pred = 0, 0
        # report = "rxn_id\tprot_id\tsubunit_id\tnumAddedGenes\tAdded\tnumRemovedGenes\tremoved\n"
        report = f"rxn_id\tprot_id\tsubunit_id\tathaFeatures\t{self.species}Features\n"
        # report = f"rxn_id\tathaFeatures\t{self.species}Features\n"

        print(bcolors.PROG+f"    Reading XML model: {self.model_name_json.replace('json', 'xml')}."+bcolors.ENDC)
        co_path = os.path.join(self.models_folder, self.model_name_json.replace('json', 'xml'))
        if os.path.exists(co_path):
            self.co_model = read_sbml_model(co_path)
        else:
            from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
            KBOF = KBaseObjectFactory()
            model_path = os.path.join(self.models_folder, self.atha_model_json)
            media_path = os.path.join(self.models_folder, "PlantAutotrophicMedia.json")

            self.co_model = KBOF.build_object_from_file(model_path, "KBaseFBA.FBAModel")
            co_media = KBOF.build_object_from_file(media_path, "KBaseBiochem.Media")
            self.co_model.medium = co_media
            print(self.co_model.optimize())
            # print(abc)

        # Keep a copy of the old genes to be removed from the model
        old_genes = self.co_model.genes.copy()
        old_genes_ids = set(self.co_model.genes._dict.keys())

        print(bcolors.SUBRESULT+f"Replacing model genes with orthologs"+bcolors.ENDC)
        print(bcolors.PROG+f"   processing reactions"+bcolors.ENDC)

        old_spc_features = set(self.spc_json_model.modelfeatures_dict.keys())
        self.spc_json_model.modelfeatures_dict = dict()
        num_moreone = 0
        # print(self.atha_json_model.modelreactions_dict.items())
        for rxn_id, atha_reaction in self.atha_json_model.modelreactions_dict.items():
            # print(rxn_id)
            if "fd" in rxn_id:
                print(atha_reaction)
                # print(abc)
            rxn_atha_ftrs, rxn_spc_ftrs_new, rxn_spc_ftrs_old = set(), set(), set()

            spc_reaction = self.spc_json_model.modelreactions_dict[rxn_id]
            # # Add reaction ID to string
            # rxn_summary_str += rxn_id+"\t"
            atha_ftrs, spc_ftrs = set(), set()

            # remove
            if len(atha_reaction.modelReactionProteins) > 1:
                num_moreone += 1
                print("----> greater than one: ", rxn_id, " ", num_moreone, " * ", len(atha_reaction.modelReactionProteins))

            mrp_list = [None] * len(atha_reaction.modelReactionProteins)
            ftrs = False
            for prot_index in range(0, len(atha_reaction.modelReactionProteins)):
                atha_mrp =atha_reaction.modelReactionProteins[prot_index]
                spc_mrp = atha_reaction.modelReactionProteins[prot_index].deepcopy()
                for subu_index in range(0, len(atha_mrp.modelReactionProteinSubunits)):
                    atha_mrps = atha_mrp.modelReactionProteinSubunits[subu_index]
                    spc_mrps = spc_mrp.modelReactionProteinSubunits[subu_index].deepcopy()
                    spc_mrps.feature_refs = list()
                    new_ftr_list = set()

                    if atha_mrps.feature_refs:
                        ftrs = True

                    for feature in atha_mrps.feature_refs:

                        f_id = feature.gene_id #feature.split("/")[-1]
                        # skipping organellar genes
                        if any(pref in f_id.upper() for pref in ['ATC', 'ATM']):
                            if f_id in self.ortho_mapping_dict:
                                print(rxn_id, "-- ", f_id, ": ", self.ortho_mapping_dict[f_id])

                        try:
                            f_id_orthologs = self.ortho_mapping_dict[f_id]
                            for ortho in f_id_orthologs:
                                newGenes.add(ortho)
                                new_ftr_list.add(feature.feature_ref.replace(f_id, ortho))
                        except KeyError:
                            notfound.add(f_id)
                    spc_mrps.update_feature_refs(list(new_ftr_list),
                                            self.spc_json_model.modelfeatures_dict)

                    spc_mrp.modelReactionProteinSubunits[subu_index] = spc_mrps

                    # Stats and logging
                    allModelGenes.update([ref.split('/')[-1] for ref in new_ftr_list])

                    try:
                        spc_old_ftrList_ids = {ftr.gene_id for ftr in \
                                        spc_reaction.modelReactionProteins[prot_index].\
                                        modelReactionProteinSubunits[subu_index].feature_refs}
                        new_ftr_list_ids = {ftr.split('/')[-1] for ftr in new_ftr_list}


                        rxn_spc_ftrs_new.update(new_ftr_list_ids)
                        rxn_spc_ftrs_old.update(spc_old_ftrList_ids)

                        new = list(new_ftr_list_ids - spc_old_ftrList_ids)
                        rmvd = list(spc_old_ftrList_ids - new_ftr_list_ids)

                    except IndexError:
                        new = [ref.split('/')[-1] for ref in new_ftr_list]
                        rmvd = []

                    # if (not new_ftr_list_ids) and spc_old_ftrList_ids:
                    #     no_pred_in_new.add(rxn_id)
                    # if (not spc_old_ftrList_ids) and new_ftr_list_ids:
                    #     no_pred_in_old.add(rxn_id)


                    allNew.update(new)
                    allRemoved.update(rmvd)
                    # report += f"{rxn_id}\t{prot_index+1}\t{subu_index+1}\t{len(new)}\t{new}\t{len(rmvd)}\t{rmvd}\n"
                    # Add atha features to string
                    atha_ftrs.update({ftr.gene_id for ftr in atha_mrps.feature_refs})
                    spc_ftrs.update({ftr.gene_id for ftr in spc_mrps.feature_refs})
                    rxn_atha_ftrs.update({ftr.gene_id for ftr in atha_mrps.feature_refs})
                    # atha_ftrs = {ftr.gene_id for ftr in atha_mrps.feature_refs}
                    # spc_ftrs  = {ftr.gene_id for ftr in spc_mrps.feature_refs}
                    # atha_ftrs = atha_ftrs if atha_ftrs else 'None'
                    # spc_ftrs  = spc_ftrs if spc_ftrs else 'None'
                    #
                    # report += f"{rxn_id}\t{prot_index+1}\t{subu_index+1}\t{atha_ftrs}\t{spc_ftrs}\n"

                mrp_list[prot_index] = spc_mrp

            if (not rxn_spc_ftrs_new) and rxn_spc_ftrs_old:
                no_pred_in_new.add(rxn_id)
            if (not rxn_spc_ftrs_old) and rxn_spc_ftrs_new:
                no_pred_in_old.add(rxn_id)
            if (not rxn_spc_ftrs_old) and (not rxn_spc_ftrs_new) and rxn_atha_ftrs:
                no_pred_at_all.add(rxn_id)
            # atha_ftrs = atha_ftrs if atha_ftrs else 'None'
            # spc_ftrs  = spc_ftrs if spc_ftrs else 'None'
            # report += f"{rxn_id}\t{atha_ftrs}\t{spc_ftrs}\n"
            spc_reaction.modelReactionProteins = mrp_list

            if ftrs:
                rxns_atha_pred += 1
            # Generate the rule for the XML model
            reaction_dict = json.loads(spc_reaction.__str__())
            rule = self.generate_gpr(reaction_dict['modelReactionProteins'], "modelReactionProteins")
            try:
                co_reaction = self.co_model.reactions.get_by_id(rxn_id)
                co_reaction.gene_reaction_rule = rule
            except KeyError:
                print(f" ----> Reaction {rxn_id} not found in co_model")

            if rule != '':
                spc_pred += 1
                # print(f"{rxn_id} {rule}")


        with open('changed_reactions.tsv','w') as fh:
            fh.write(report)
        # remove old genes from model
        print(bcolors.PROG+f"   cleaning up sbml model"+bcolors.ENDC)
        new_genes = set()
        for gene in allModelGenes:
            new_genes.add(self.co_model.genes.get_by_id(gene))

        remove_genes = set(old_genes) - set(new_genes)
        print(bcolors.PROG+f"   removing {len(remove_genes)} (model had {len(old_genes_ids)} now there are {len(new_genes)})"+bcolors.ENDC)
        co.manipulation.remove_genes(model=self.co_model, gene_list=remove_genes, remove_reactions=False)

        print(f"There are {rxns_atha_pred} A. Thaliana reactions with predicted genes and {spc_pred} in {self.species}")
        if vebose:
            print(bcolors.WARNING+f"No orthologs found for {len(notfound)} genes"+bcolors.ENDC)

            print(bcolors.WARNING+f"New model genes {len(allNew-old_genes_ids)}"+bcolors.ENDC)
            # print(bcolors.PROG+f"---> {allNew} ")
            print(bcolors.WARNING+f"Removed model genes {len(allRemoved-allModelGenes)}"+bcolors.ENDC)
            # print(bcolors.PROG+f"---> {allRemoved-allModelGenes} ")
            print(bcolors.WARNING+f"{len(no_pred_in_new)} Reactions w/out predictions in new model, but w/ predictions in old model. "+bcolors.ENDC)
            print(bcolors.PROG+f"---> {no_pred_in_new} ")

            print(bcolors.WARNING+f"{len(no_pred_in_old)} reactions with genes in the new model that did not have any in the old model "+bcolors.ENDC)
            print(f"No prediction in old or new model: {len(no_pred_at_all)}")
            print(bcolors.PROG+f"---> {no_pred_at_all} ")

            self.allFeatures = set(self.co_model.genes._dict.keys())
            print(f"---> coPy genes {len(self.allFeatures)}")

    def saveModel2File(self, sbicolor5=False):
        model_name_json = self.model_name_json
        # model_name_json = self.model_name_json.replace(".json", \
                                # f"_new_orthologs_{self.version}.json")
        if sbicolor5:
            model_name_json = self.model_name_json.replace("3.1.1", "5.1")
            self.spc_json_model.name = self.spc_json_model.name.replace("3.1.1", "5.1")
            self.spc_json_model.id   = self.spc_json_model.id.replace("3.1.1", "5.1")

        if sbicolor5 and ("5.1" not in model_name_json):
            model_name_json = model_name_json.replace('.json', '5.1.json')

        print(bcolors.PROG+f"   Saving new models to {model_name_json}/.xml"+bcolors.ENDC)
        json_output = json.dumps(self.spc_json_model.as_dict(), indent=4)

        if not os.path.exists(os.path.join(self.models_folder, f"ortho_{self.version}_models")):
            os.makedirs(os.path.join(self.models_folder, f"ortho_{self.version}_models"))

        output_json_file = os.path.join(self.models_folder, f"ortho_{self.version}_models", "new_"+model_name_json)
        with open(output_json_file, 'w') as outfile:
            outfile.write(json_output)

        from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
        KBOF = KBaseObjectFactory()
        # model_path = os.path.join(self.models_folder, self.atha_model_json)
        media_path = os.path.join(self.models_folder, "PlantAutotrophicMedia.json")

        self.co_model = KBOF.build_object_from_file(output_json_file, "KBaseFBA.FBAModel")
        co_media = KBOF.build_object_from_file(media_path, "KBaseBiochem.Media")
        self.co_model.medium = co_media
        # print(self.co_model.optimize())
        write_sbml_model(self.co_model, output_json_file.replace(".json", ".xml"))
        solution = self.co_model.optimize()
        print('----------->> ', solution)
        return output_json_file

    def reactionsWithGenes(self):
        total = 0

        print(bcolors.PROG+f"    Reading XML model: {self.model_name_json.replace('json', 'xml')}."+bcolors.ENDC)
        # self.co_model = read_sbml_model(os.path.join(self.models_folder, self.model_name_json.replace('json', 'xml'))) #atha_model_json
        self.co_model = read_sbml_model("/Users/sea/Projects/QPSI_project/QPSI_Modeling/data/metabolic_models/fullmodels_media/sbicolor_3.1.1_modelnew_orthologs.xml")

        for rxn in self.co_model.reactions:
            if rxn.gene_reaction_rule != '':
                total += 1
        print(f"{self.species} {total}/{len(self.co_model.reactions)}")

    def modelStats(self, output_json_file):
        #### READ ROLES -----------
        PS_url  = "https://raw.githubusercontent.com/ModelSEED/PlantSEED/"
        PS_tag  = "8cf60046e4af68912f7a7d3eeff16880a07f56bd"
        PS_json = "/Data/PlantSEED_v3/PlantSEED_Roles.json"

        print(bcolors.SUBRESULT+"PS Roles reading ...")
        PS_json_data = json.load(urlopen(PS_url+PS_tag+PS_json))
        roles_dict = dict()
        for item in PS_json_data:
            # Create a new role entry for the consolidated model
            role_dict = dict()
            role_dict["role"]        = item["role"]
            role_dict["subsystems"]  = item["subsystems"]
            role_dict["reactions"]   = list()
            roles_dict[item["role"]] = role_dict

        # READ MODEL -----------
        metModel = None
        with open(output_json_file, 'r') as f:
            data=f.read()
            json_model = json.loads(data)
            metModel = mc.Model().fromJSON(json_model, roles_dict)
            if metModel == None:
                return

        # PROCESS REACTIONS
        subsys_dict = dict() # reactions - wGenes - wtGenes
        transporters = 0
        all_roles = set()
        for rxnId, rxn in metModel.modelreactions_dict.items():
            if rxn.isTransporter() != 'NotTransport':
                transporters += 1
                continue

            all_roles.update(set(rxn.role_properties().keys()))

            for subsys in rxn.subsystems:
                try:
                    sub_dict = subsys_dict[subsys]
                except KeyError:
                    sub_dict = {'reactions':set(), 'wGenes':0, 'wtGenes':0}

                if rxnId in sub_dict['reactions']:
                    continue
                else:
                    sub_dict['reactions'].add(rxnId)
                    if rxn.genes:
                        sub_dict['wGenes'] = sub_dict['wGenes']+1
                    else:
                        sub_dict['wtGenes'] = sub_dict['wtGenes']+1

                subsys_dict[subsys] = sub_dict

        notProcessed = list(subsys_dict.keys())
        pathways_class_dict = dict()
        for item in PS_json_data:
            all_classes = item["classes"].keys()
            for cls in all_classes:
                for pathway in item["classes"][cls]:
                    if pathway in notProcessed:
                        pathways_class_dict[pathway] = cls
                        notProcessed.remove(pathway)
                    if not notProcessed: break

        dict_str = "subsystem\tclass\tnumReactions\twGenes\twtGenes\n"
        for subsys, sdict in subsys_dict.items():
            cls = pathways_class_dict[subsys]
            subsys = subsys.replace("_", " ").replace(" in plants", "")
            dict_str += "{}\t{}\t{}\t{}\t{}\n".format(subsys, cls, \
                    len(sdict['reactions']), sdict['wGenes'], sdict['wtGenes'])


        with open(f'subsystemStats_{self.species}.tsv','w') as fh:
            fh.write(dict_str)

        print(f"There are {len(metModel.modelreactions_dict)} reactions.")
        print(f" --- {transporters} transporters.")
        print(f"There are {len(metModel.modelcompounds_dict)} metabolites.")
        print(f"There are {len(metModel.modelfeatures_dict)} model genes.")
        print(f"There are {len(all_roles)} functional roles.")
        print(f"There are {len(subsys_dict)} metabolic pathways.")

        # GENE STATS
        OG_set = set()
        oos, pah, pao = 0, 0, 0
        for gene, val in self.ortho_dict.items():
            if gene not in metModel.modelfeatures_dict: continue
            OG_set.add(val[0])
            if val[1] == 'O': oos += 1
            if val[1] == 'PAH': pah += 1
            if val[1] == 'PAO': pao += 1

        # print(self.ortho_mapping_dict)
        print(f"There are {oos} O's")
        print(f"There are {pah} PAH's")
        print(f"There are {pao} PAO's")
        print(f"There are {len(OG_set)} orthologous groups")


if __name__ == '__main__':
    mode = 3
    ortho_version = 'apr10'

    # Load the model (replace 'your_model.xml' with your actual file path)
    model = read_sbml_model("/Users/selalaoui/Projects/QPSI_project/paper_repo/RNASeq_Enzyme_Abundance/data/metabolic_models/fullmodels_media/Athaliana_Thylakoid_Reconstruction_ComplexFix_070224.xml")

    # Print the number of genes in the model
    print(f"Number of genes in atha model: {len(model.genes)}")

    model = read_sbml_model("/Users/selalaoui/Projects/QPSI_project/paper_repo/RNASeq_Enzyme_Abundance/data/metabolic_models/fullmodels_media/ptrich_4.1_Thylakoid_Reconstruction_ComplexFix_070224.xml")

    # Print the number of genes in the model
    print(f"Number of genes in poplar model: {len(model.genes)}")


    model = read_sbml_model("/Users/selalaoui/Projects/QPSI_project/paper_repo/RNASeq_Enzyme_Abundance/data/metabolic_models/plastidial_models/Athaliana_plastid_Thylakoid_Reconstruction_ComplexFix_070224.xml")

    # Print the number of genes in the model
    print(f"Number of genes in atha plastid model: {len(model.genes)}")

    model = read_sbml_model("/Users/selalaoui/Projects/QPSI_project/paper_repo/RNASeq_Enzyme_Abundance/data/metabolic_models/plastidial_models/ptrich_4.1_plastid_Thylakoid_Reconstruction_ComplexFix_RevFix3_250512.xml")

    # Print the number of genes in the model
    print(f"Number of genes in poplar plastid model: {len(model.genes)}")
    
    print(abc)

    spc_version = '3.1.1' # Sorghum: '3.1.1' Poplar: '4.1'
    spc = "sbicolor" # ['ptrich', 'sbicolor']

    spc_version = '4.1' # Sorghum: '3.1.1', 5.1 Poplar: '4.1'
    spc = "ptrich"

    sbicolor5=('5' in spc_version)

    mdpHelperObj = ModelDataProcessingHelper(species=spc, ortho_version=ortho_version)#
    mdpHelperObj.readModelJSON()
    # mdpHelperObj.reactionsWithGenes()

    mdpHelperObj.readOrthologs(mode=mode, version=spc_version)
    mdpHelperObj.replaceOrthologsInModel()

    output_json_file = mdpHelperObj.saveModel2File(sbicolor5=sbicolor5)
    mdpHelperObj.modelStats(output_json_file)
