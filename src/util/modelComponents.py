import json
import re
import copy
import os
import sys
from typing import TYPE_CHECKING, Dict, List, TypeVar, Generic
T = TypeVar('T')

## ----------------------------------------------------------------------------------------
# Collection of classes representing all the elements of the metabolic model.
# Each class provides helper methods to ease the processing of the model and extraction of
#   model components.
# The classes are used for plastidial model extraction and for reaction expression scores
#   computation
## ----------------------------------------------------------------------------------------

class Util(Generic[T]):
    def addReactionReagents(self, rxn_reagents, new_reagents_list:List[T]):
        if new_reagents_list == None:
            raise ValueError(bcolors.FAIL+"  No reaction reagents provided."+bcolors.ENDC)

        if not isinstance(new_reagents_list, list) or (not all(isinstance(r, self.__orig_class__.__args__[0]) for r in new_reagents_list)):
            raise ValueError(bcolors.FAIL+"  Wrong Type -- reaction reagents not added."+bcolors.ENDC)

        if rxn_reagents == None:
            rxn_reagents = new_reagents_list
        else:
            rxn_reagents.extend(new_reagents_list)

        return rxn_reagents

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

class Compartment:
    def __init__(self, cpmt_letterId:str, cpmt_index:int, label:str=None, pH:int=0, potential:int=0):
        self.compartmentIndex =  cpmt_index
        self.compartment_ref =  "~/template/compartments/id/"+cpmt_letterId
        self.id =  cpmt_letterId+str(cpmt_index)
        self.label =  label
        self.pH =  pH
        self.potential = potential

    @classmethod
    def fromJSON(cls, cpmt_dict: Dict) -> None:
        """creat `cls` from JSON dict """
        letter_id = cpmt_dict["compartment_ref"].split("/")[-1]
        return cls(letter_id, cpmt_dict["compartmentIndex"], cpmt_dict["label"],
                    cpmt_dict["pH"], cpmt_dict["potential"])

    def deepcopy(self):
        return copy.deepcopy(self)

    def __str__(self) -> str:
        return json.dumps(self.__dict__, sort_keys=False, indent=4)

class Compound:
    def __init__(self, id:str, name:str, formula:str, charge:float=0.0):
        if "_" not in id:
            raise ValueError("Invalid compound ID")

        self.charge =  charge
        self.compound_ref =  "~/template/compounds/id/"+id.split('_')[0]
        self.formula =  formula
        self.id =  id
        self.modelcompartment_ref =  "~/modelcompartments/id/"+id.split('_')[-1]
        self.name =  name

    @classmethod
    def fromJSON(cls, cpd_dict: Dict) -> None:
        """creat `cls` from JSON dict """
        return cls(cpd_dict["id"], cpd_dict["name"], cpd_dict["formula"],
                    cpd_dict["charge"])

    @property
    def compartment(self):
        return self.id.split('_')[-1]

    def update_cpd_compartment(self, new_cprmt:str) -> None:
        # update ID
        self.id = re.sub(r"_[a-z]0", new_cprmt, self.id)
        # update comparment ref
        self.modelcompartment_ref = "~/modelcompartments/id/"+new_cprmt.replace('_', '')
        # update compound ref
        self.compound_ref = "~/template/reactions/id/"+self.id.split('_')[0]

    def deepcopy(self):
        return copy.deepcopy(self)

    def __str__(self) -> str:
        return json.dumps(self.__dict__, sort_keys=False, indent=4)

class Reagent:
    def __init__(self, coeff:int = 0, cpd_id:str = ""):
        self.coefficient = coeff
        if "modelcompounds" in cpd_id:
            self.modelcompound_ref = cpd_id
        else:
            self.modelcompound_ref = "~/modelcompounds/id/"+cpd_id

    @property
    def id(self) -> str:
        return self.modelcompound_ref.split('/')[-1]

    @property
    def compartment(self) -> str:
        return self.modelcompound_ref.split('_')[-1]

    def update_reagent_compartment(self, new_cprmt:str) -> None:
        self.modelcompound_ref = self.modelcompound_ref.replace("_"+self.compartment, new_cprmt)

    @classmethod
    def fromJSON(cls, rgt_dict: Dict) -> None:
        """creat `cls` from JSON dict """
        return cls(rgt_dict["coefficient"], rgt_dict["modelcompound_ref"])

    def toDict(self):
        return self.__dict__

    def __str__(self) -> str:
        return json.dumps(self.__dict__, sort_keys=False, indent=4)

class Feature:
    def __init__(self, feature_ref:str, weight=-1.0, metalBinding:set=set(), pfam:set=set()):
        self.feature_ref  = feature_ref
        self.weight       = weight
        self.binds        = metalBinding
        self.pfam         = pfam

    @property
    def gene_id(self):
        return self.feature_ref.split("/")[-1]

    @property
    def binds(self):
        return self._binds

    @binds.setter
    def binds(self, binds):
        if not hasattr(self, "_binds"):
            self._binds= set()

        # print(f"-- Setting binds {binds}")
        if isinstance(binds, str):
            self._binds.add(binds)
        else:
            self._binds.update(binds)


    # @property
    # def pfam(self):
    #     return self.pfam
    #
    # @pfam.setter
    # def pfam(self, pfam:str):
    #     self._pfam = pfam
    #
    # @property
    # def weight(self):
    #     return self.weight
    #
    # @weight.setter
    # def weight(self, weight:float):
    #     self._weight = weight

    def updateBinding(self, newMatel):
        if isinstance(newMatel, str):
            self.metalBinding.add(newMatel)
        else:
            self.metalBinding.update(newMatel)

    def deepcopy(self):
        return copy.deepcopy(self)


    def __str__(self) -> str:
        return f"{self.feature_ref}"

class ProteinSubUnit:
    def __init__(self, feature_refs:List[Feature]=None, note:str='Features uncharacterized and unannotated', optionalSubunit:int=0, role:str='Unknown', triggering:int=0, subsystems:set=set()):
        self.feature_refs    = feature_refs
        self.note            = note
        self.optionalSubunit = optionalSubunit
        self.role            = role
        self.triggering      = triggering
        self.subsystems      = subsystems

    def deepcopy(self) -> "ProteinSubUnit":
        return copy.deepcopy(self)

    def toDict(self):
        temp_dict = self.__dict__.copy()
        temp_dict['feature_refs'] = [ftr.__str__()
                                for ftr in temp_dict['feature_refs']]
        temp_dict.pop('subsystems')
        return temp_dict

    def update_feature_refs(self, new_feature_refs:List[str]=None, model_ftr_dict=dict()):
        if not new_feature_refs:
            return

        self.feature_refs = list()
        for ftr_ref in new_feature_refs:
            ftr = Feature(ftr_ref)
            if model_ftr_dict:
                ftr_ref_id = ftr_ref.split('/')[-1]
                if ftr_ref_id in model_ftr_dict:
                    ftr = model_ftr_dict[ftr_ref_id]
                else:
                    model_ftr_dict[ftr_ref_id] = ftr
            self.feature_refs.append(ftr)


    def __str__(self) -> str:
        return json.dumps(self.toDict(), sort_keys=False, indent=4)

class ReactionProtein:
    def __init__(self, complex_ref:str="", modelReactionProteinSubunits:List[ProteinSubUnit]=None, note:str="None"):
        self.complex_ref = complex_ref
        self.modelReactionProteinSubunits = modelReactionProteinSubunits
        self.note = note

    @classmethod
    def fromJSON(cls, rp_dict: Dict, ftr_dict:Dict=dict(), roles_dict:Dict=dict()):
        """creat `cls` from JSON dict """
        psu_list = list()
        if "modelReactionProteinSubunits" in rp_dict:
            for rpsu in rp_dict["modelReactionProteinSubunits"]:
                # Create the list of features objects
                features = list()
                for ftr_ref in rpsu["feature_refs"]:
                    gene_id = ftr_ref.split('/')[-1]
                    try:
                        ftr = ftr_dict[gene_id]
                    except KeyError:
                        ftr = Feature(ftr_ref)
                        ftr_dict[gene_id] = ftr
                    features.append(ftr)
                # Retrieve the subsystems
                subsystems = set()
                if roles_dict:
                    if rpsu["role"] in roles_dict:
                        subsystems = roles_dict[rpsu["role"]]['subsystems']

                # Create new subunit
                psu_list.append(ProteinSubUnit(features, rpsu["note"],
                rpsu["optionalSubunit"], rpsu["role"], rpsu["triggering"], subsystems))

        if len(psu_list) == 0:
            psu = ProteinSubUnit([], note="Features uncharacterized and unannotated")
            psu_list = [psu]

        complex_ref = rp_dict["complex_ref"] if "complex_ref" in rp_dict else ""
        note = rp_dict["note"] if "note" in rp_dict else ""

        return cls(complex_ref, psu_list, note)

    @property
    def reactionProtein(self):
        return self.__str__

    @reactionProtein.setter
    def reactionProtein(self, reactionProteins_dict: Dict, ftr_dict:Dict=dict(), roles_dict:Dict=dict()) -> None:
        temp_prot = self.fromJSON(reactionProteins_dict, ftr_dict, roles_dict)
        self.complex_ref                  = temp_prot.complex_ref
        self.modelReactionProteinSubunits = temp_prot.modelReactionProteinSubunits
        self.note                         = temp_prot.note

    def toDict(self):
        temp_dict = self.__dict__.copy()
        temp_dict['modelReactionProteinSubunits'] = [psu.toDict()
                                for psu in temp_dict['modelReactionProteinSubunits']]
        return temp_dict

    def deepcopy(self) -> "ReactionProtein":
        return copy.deepcopy(self)

    def __str__(self) -> str:
        return json.dumps(self.toDict(), sort_keys=False, indent=4)
        # return json.dumps(self.__dict__, default=self.to_serializable, sort_keys=False, indent=4)

class Reaction:
    def __init__(self, id:str, direction:str, modelReactionReagents:List[Reagent],
     modelReactionProteins:List[ReactionProtein]=None, name:str="",
     probability:int=0, protons:int=0):
        if "_" not in id:
            raise ValueError("Invalid reaction ID")

        self.direction             = direction
        self.id                    = id
        self.modelReactionProteins = list()
        self.modelReactionReagents = list()
        self.modelcompartment_ref  = "~/modelcompartments/id/"+id.split('_')[-1]
        self.name                  = name
        self.probability           = probability
        self.protons               = protons
        self.reaction_ref          = "~/template/reactions/id/"+id[:-1]

        # protein_list = Util[Reagent]().addReactionReagents(self.modelReactionProteins, modelReactionProteins)
        self.addReactionProteins(modelReactionProteins)
        self.modelReactionReagents = Util[Reagent]().addReactionReagents(self.modelReactionReagents, modelReactionReagents)

    def addReactionProteins(self, protein_list:List[ReactionProtein]):
        if isinstance(protein_list, list) and (all(isinstance(p, ReactionProtein) for p in protein_list)):
            if self.modelReactionProteins == None:
                self.modelReactionProteins = protein_list
            else:
                self.modelReactionProteins.extend(protein_list)
        else:
            print(bcolors.PROG+"\t    {} Wrong Type or protein list empty -- Adding empty list instead.".format(self.id)+bcolors.ENDC)
            psu = ProteinSubUnit([], note="Features uncharacterized and unannotated")

            self.modelReactionProteins = [ReactionProtein(modelReactionProteinSubunits=[psu], note="universal")]

    @classmethod
    def fromJSON(cls, rxn_dict: Dict, ftr_dict:Dict=dict(), roles_dict:Dict=dict()):
        """creat `cls` from JSON dict """
        r_list = [Reagent.fromJSON(e) for e in rxn_dict["modelReactionReagents"]]
        p_list = [ReactionProtein.fromJSON(e, ftr_dict, roles_dict) for e in rxn_dict["modelReactionProteins"]]

        return cls(rxn_dict["id"], rxn_dict["direction"], r_list, p_list,
        rxn_dict["name"], rxn_dict["probability"], rxn_dict["protons"])

    @property
    def idOnly(self):
        return self.id.split('_')[0]

    @property
    def compartment(self):
        return self.id.split('_')[-1]

    @property
    def compartments(self) -> set:
        return {r.compartment for r in self.modelReactionReagents}

    @property
    def reagents(self) -> set:
        return {r.id for r in self.modelReactionReagents}

    @property
    def reactants(self) -> set:
        # return all if rxn is reversible
        if self.direction == '=': return self.reagents
        # return reagents with negative coeff
        return {r.id for r in self.modelReactionReagents if r.coefficient < 0}

    @property
    def products(self) -> set:
        # return all if rxn is reversible
        if self.direction == '=': return self.reagents
        # return reagents with positive coeff
        return {r.id for r in self.modelReactionReagents if r.coefficient >= 0}

    @property
    def genes(self) -> set:
        genes = set()
        for mrp in self.modelReactionProteins:
            for mpsu in mrp.modelReactionProteinSubunits:
                for ftr in mpsu.feature_refs:
                    genes.add(ftr.gene_id)

        return genes

    @property
    def gpr(self):
        gpr = list()
        for mrp in self.modelReactionProteins:
            mpsu_list = list()
            for mpsu in mrp.modelReactionProteinSubunits:
                mpsu_list.append([ftr.gene_id for ftr in mpsu.feature_refs])
            gpr.append(mpsu_list)

        return gpr

    @property
    def subsystems(self):
        subsystems = set()
        for mrp in self.modelReactionProteins:
            for mpsu in mrp.modelReactionProteinSubunits:
                subsystems.update(mpsu.subsystems)

        return subsystems

    @property
    def binds(self):
        binds = set()
        for mrp in self.modelReactionProteins:
            for mpsu in mrp.modelReactionProteinSubunits:
                for ftr in mpsu.feature_refs:
                    binds.update(ftr.binds)

        return binds

    def role_properties(self, metal='none'):
        roles_dict     = dict()

        for mrp in self.modelReactionProteins:
            for mpsu in mrp.modelReactionProteinSubunits:
                ftr_list      = list()
                binding_genes = set()
                pfam_set      = set()
                bind          = "none"
                for ftr in mpsu.feature_refs:
                    ftr_list.append(ftr.gene_id)
                    pfam_set.update(ftr.pfam)
                    if metal in ftr.binds:
                        bind = metal
                        binding_genes.add(ftr.gene_id)

                roles_dict[mpsu.role] = {"features": ftr_list,
                                        "predicted_genes": binding_genes,
                                        "PFAM_domains": pfam_set,
                                        "bind": bind}

        return roles_dict

    def update_reaction_compartment(self, new_cprmt:str, is_trans:str) -> None:
        curr_cprmt = self.compartment
        # update ID
        self.id = re.sub(r"_[a-z]0", new_cprmt, self.id)
        # update comparment ref
        self.modelcompartment_ref = "~/modelcompartments/id/"+new_cprmt.replace('_', '')
        # update reaction ref
        self.reaction_ref = "~/template/reactions/id/"+self.id[:-1]
        # update reagents compartments
        for r in self.modelReactionReagents:
            if is_trans!='NotTransport' and r.compartment == curr_cprmt:
                r.update_reagent_compartment(new_cprmt)
            else:
                r.update_reagent_compartment(new_cprmt)

    def deepcopy(self) -> "Reaction":
        return copy.deepcopy(self)

    def isTransporter(self) -> str:
        compartment_set = self.compartments
        if len(compartment_set) == 1:
            return 'NotTransport'
        elif 'e0' in compartment_set:
            return 'MediumTransport'
        else:
            return 'OtherTransport'

    @staticmethod
    def to_serializable(value):
        if hasattr(value, '__dict__'):
            return value.__dict__
        return value

    def toDict(self):
        temp_dict = self.__dict__.copy()
        temp_dict['modelReactionProteins'] = [mrp.toDict()
                            for mrp in temp_dict['modelReactionProteins']]
        temp_dict['modelReactionReagents'] = [rgt.toDict()
                            for rgt in temp_dict['modelReactionReagents']]
        return temp_dict

    def __str__(self) -> str:
        return json.dumps(self.toDict(), sort_keys=False, indent=4)

class Biomass:
    def __init__(self, id:str, name:str, biomasscompounds:List[Reagent], cellwall:int=0,
    cofactor:int=0, dna:int=0, energy:int=0, lipid:int=0, other:int=0,
    protein:int=0, rna:int=0):
        self.cellwall = cellwall
        self.cofactor = cofactor
        self.dna = dna
        self.energy = energy
        self.id = id
        self.lipid = lipid
        self.name = name
        self.other = other
        self.protein = protein
        self.rna = rna

        self.biomasscompounds = None
        self.biomasscompounds = Util[Reagent]().addReactionReagents(self.biomasscompounds, biomasscompounds)

    @classmethod
    def fromJSON(cls, bio_dict: Dict):
        """creat `cls` from JSON dict """
        cpd_list = [Reagent.fromJSON(e) for e in bio_dict["biomasscompounds"]]

        return cls(bio_dict["id"], bio_dict["name"], cpd_list,
        bio_dict["cellwall"], bio_dict["cofactor"], bio_dict["dna"],
        bio_dict["energy"], bio_dict["lipid"], bio_dict["other"],
        bio_dict["protein"], bio_dict["rna"])

    @staticmethod
    def to_serializable(value):
        if hasattr(value, '__dict__'):
            return value.__dict__
        return value

    def __str__(self) -> str:
        return json.dumps(self.__dict__, default=self.to_serializable, sort_keys=False, indent=4)

class Model:
    def __init__(self, name:str=None, id:str=None, modelReactions:dict=dict(),
                biomasses:dict=dict(), modelCompartments:dict=dict(),
                modelCompounds:dict=dict(), source:str=None, source_id:str=None,
                gapfillings:list=None, gapgens:list=None, template_ref:str=None,
                genome_ref:str=None, type:str=None, modelfeatures:dict=dict()):
        self.name                   = name
        self.id                     = id
        self.source                 = source
        self.source_id              = source_id
        self.gapfillings            = gapfillings
        self.gapgens                = gapgens
        self.template_ref           = template_ref
        self.genome_ref             = genome_ref
        self.type                   = type
        self.modelreactions_dict    = modelReactions
        self.biomasses_dict         = biomasses
        self.modelcompartments_dict = modelCompartments
        self.modelcompounds_dict    = modelCompounds
        self.modelfeatures_dict     = modelfeatures

    @classmethod
    def fromJSON(cls, model: Dict, roles_dict:Dict=dict()):
        """creat `cls` from JSON dict """

        # rxn_dict, ftr_dict = self.reactionFeaturesFromJSON(model["modelreactions"])
        ftr_dict = dict()
        rxn_dict = {r['id']:Reaction.fromJSON(r, ftr_dict, roles_dict) for r in model["modelreactions"]}
        # print(len(ftr_dict))
        # print(ftr_dict)
        cpmt_dict = {c['id']:Compartment.fromJSON(c) for c in model["modelcompartments"]}
        cpd_dict = {c['id']:Compound.fromJSON(c) for c in model["modelcompounds"]}
        bio_dict = {b['id']:Biomass.fromJSON(b) for b in model["biomasses"]}

        return cls(model["name"], model["id"], rxn_dict, bio_dict, cpmt_dict,
                cpd_dict, model["source"], model["source_id"], model["gapfillings"],
                model["gapgens"], model["template_ref"], model["genome_ref"], model["type"], ftr_dict)

    def set_model_attributes(self, name:str=None, id:str=None, source:str=None,
                source_id:str=None, gapfillings:list=None, gapgens:list=None,
                template_ref:str=None, genome_ref:str=None, type:str=None) -> None:
        self.name = name
        self.id = id
        self.source = source
        self.source_id = source_id
        self.gapfillings = gapfillings
        self.gapgens = gapgens
        self.template_ref = template_ref
        self.genome_ref = genome_ref
        self.type = type


    def add_metabolite_toCprmt(self, m_id, curr_cprmt, new_cprmt):
        if m_id not in self.modelcompounds_dict:
            curr_id = m_id.replace(new_cprmt, curr_cprmt)
            try:
                new_cpd = self.modelcompounds_dict[curr_id].deepcopy()
                new_cpd.update_cpd_compartment(new_cprmt)
                self.modelcompounds_dict[m_id] = new_cpd
            except KeyError:
                print(bcolors.PROG+"Skipping {}: metabolite not found in model compounds.".format(curr_id)+bcolors.ENDC)

    def add_new_metabolite(self, metabolite) -> None:
        if isinstance(metabolite, Compound):
            self.modelcompounds_dict[metabolite.id] = metabolite
        elif isinstance(metabolite, list) and any(isinstance(m, Compound) for m in metabolite):
            for m in metabolite:
                if isinstance(m, Compound): self.modelcompounds_dict[m.id] = m
        else:
            print(bcolors.WARNING+"Trying to add metabolite(s) to model ...\n Wrong type, skipping."+bcolors.ENDC)

    def add_new_reaction(self, reaction:Reaction) -> None:
        if isinstance(reaction, Reaction):
            self.modelreactions_dict[reaction.id] = reaction
        else:
            print(bcolors.WARNING+"Trying to add reaction to model ...\n Wrong type, skipping."+bcolors.ENDC)

    def add_new_compartment(self, cprmt:Compartment) -> None:
        if isinstance(cprmt, Compartment):
            self.modelcompartments_dict[cprmt.id] = cprmt
        else:
            print(bcolors.WARNING+"Trying to add compartment to model ...\n Wrong type, skipping."+bcolors.ENDC)

    def add_new_biomass(self, bio:Biomass) -> None:
        if isinstance(bio, Biomass):
            self.biomasses_dict[bio.id] = bio
        else:
            print(bcolors.WARNING+"Trying to add biomass reaction to model ...\n Wrong type, skipping."+bcolors.ENDC)

    def compute_atoms_franctions(self, bio_id:str=None):
        if self.biomasses_dict == None or len(self.biomasses_dict) == 0:
            print(bcolors.WARNING+"\t Biomass not initialized in model yet."+bcolors.ENDC)
            return

        from BiochemPy import Compounds
        cpd_helper = Compounds()
        atom_dict = dict()

        try:
            bio = self.biomasses_dict[bio_id]
        except KeyError:
            bio = list(self.biomasses_dict.values())[0]

        for cpd in bio.biomasscompounds:
            # Get coefficient and formula
            coeff = cpd.coefficient
            cpd_formula = self.modelcompounds_dict[cpd.id].formula

            # Extract the atoms and add to add to dict
            cpd_atoms = cpd_helper.parseFormula(cpd_formula)
            for atom, stoichiometry in cpd_atoms.items():
                if atom in atom_dict:
                    atom_dict[atom]= atom_dict[atom]+ stoichiometry*coeff
                else:
                    atom_dict[atom] = stoichiometry*coeff
        return atom_dict

    @staticmethod
    def to_serializable(value):
        if hasattr(value, '__dict__'):
            return value.__dict__
        return value

    # return dict values NOT the dict itself
    def as_dict(self) -> dict:
        model_dict = dict()
        model_dict["name"] = self.name
        model_dict["id"] = self.id
        model_dict["source"] = self.source
        model_dict["source_id"] = self.source_id
        model_dict["gapfillings"] = self.gapfillings
        model_dict["gapgens"] = self.gapgens
        model_dict["template_ref"] = self.template_ref
        model_dict["genome_ref"] = self.genome_ref
        model_dict["type"] = self.type
        # Sort all other elements by ID
        model_dict["modelreactions"] = [json.loads(r.__str__()) for r in sorted(self.modelreactions_dict.values(), key=lambda rxn:rxn.id)]
        model_dict["biomasses"] = [json.loads(b.__str__()) for b in sorted(self.biomasses_dict.values(), key=lambda bio:bio.id)]
        model_dict["modelcompartments"] = [json.loads(c.__str__()) for c in sorted(self.modelcompartments_dict.values(), key=lambda cpt:cpt.id)]
        model_dict["modelcompounds"] = [json.loads(cpd.__str__()) for cpd in sorted(self.modelcompounds_dict.values(), key=lambda cpd:cpd.id)]

        return model_dict

if __name__ == '__main__':
    rd = {
        "direction": ">",
        "id": "rxn00154_e0",
        "modelReactionProteins": [
            {
                "complex_ref": "~/template/complexes/id/Cpx.63",
                "modelReactionProteinSubunits": [
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.001G198000",
                            "~/genome/features/id/Potri.003G043900",
                            "~/genome/features/id/Potri.008G027400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide acetyltransferase component of pyruvate dehydrogenase complex (EC 2.3.1.12)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.008G100800",
                            "~/genome/features/id/Potri.010G151400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide dehydrogenase of pyruvate dehydrogenase complex (EC 1.8.1.4)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.008G192500",
                            "~/genome/features/id/Potri.010G038400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Pyruvate dehydrogenase E1 component alpha subunit (EC 1.2.4.1)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.001G061400",
                            "~/genome/features/id/Potri.003G166400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Pyruvate dehydrogenase E1 component beta subunit (EC 1.2.4.1)",
                        "triggering": 1
                    }
                ],
                "note": "universal"
            }
        ],
        "modelReactionReagents": [
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00003_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00004_e0"
            },
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00010_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00011_e0"
            },
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00020_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00022_e0"
            }
        ],
        "modelcompartment_ref": "~/modelcompartments/id/e0",
        "name": "R00209",
        "probability": 0,
        "protons": 0,
        "reaction_ref": "~/template/reactions/id/rxn00154_e"
    }

    # print(rdr.compartment)
    # r1 = rdr
    # r1.update_reaction_compartment('_m0', r1.isTransporter())
    # di = json.loads(r1.__str__())
    # # print(di)

    rxn = {
        "direction": "=",
        "id": "rxn01361_c0",
        "modelReactionProteins": [
            {
                "complex_ref": "Unknown",
                "modelReactionProteinSubunits": [
                    {
                        "feature_refs": [],
                        "note": "Features uncharacterized and unannotated",
                        "optionalSubunit": 0,
                        "role": 'null',
                        "triggering": 0
                    }
                ],
                "note": "universal"
            }
        ],
        "modelReactionReagents": [
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00003_c0"
            },
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00282_c0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00004_c0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00247_c0"
            }
        ],
        "modelcompartment_ref": "~/modelcompartments/id/c0",
        "name": "R01869",
        "probability": 0,
        "protons": 0,
        "reaction_ref": "~/template/reactions/id/rxn01361_c"
    }
    rd = {
        "direction": ">",
        "id": "rxn00154_e0",
        "modelReactionProteins": [
            {
                "complex_ref": "~/template/complexes/id/Cpx.63",
                "modelReactionProteinSubunits": [
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.001G198000",
                            "~/genome/features/id/Potri.003G043900",
                            "~/genome/features/id/Potri.008G027400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide acetyltransferase component of pyruvate dehydrogenase complex (EC 2.3.1.12)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.008G100800",
                            "~/genome/features/id/Potri.010G151400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide dehydrogenase of pyruvate dehydrogenase complex (EC 1.8.1.4)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.008G192500",
                            "~/genome/features/id/Potri.010G038400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Pyruvate dehydrogenase E1 component alpha subunit (EC 1.2.4.1)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Pyruvate dehydrogenase E1 component alpha subunit (EC 1.2.4.1)",
                        "triggering": 1
                    },
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.001G061400",
                            "~/genome/features/id/Potri.003G166400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Pyruvate dehydrogenase E1 component beta subunit (EC 1.2.4.1)",
                        "triggering": 1
                    }
                ],
                "note": "universal"
            },
            {
                "complex_ref": "~/template/complexes/id/Cpx.63",
                "modelReactionProteinSubunits": [
                    {
                        "feature_refs": [],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide acetyltransferase component of pyruvate dehydrogenase complex (EC 2.3.1.12)",
                        "triggering": 1
                    }
                ],
                "note": "universal"
            },
            {
                "complex_ref": "~/template/complexes/id/Cpx.63",
                "modelReactionProteinSubunits": [
                    {
                        "feature_refs": [
                            "~/genome/features/id/Potri.001G198000",
                            "~/genome/features/id/Potri.003G043900",
                            "~/genome/features/id/Potri.008G027400"
                        ],
                        "note": "Features characterized and annotated",
                        "optionalSubunit": 0,
                        "role": "Dihydrolipoamide acetyltransferase component of pyruvate dehydrogenase complex (EC 2.3.1.12)",
                        "triggering": 1
                    }
                ],
                "note": "universal"
            }
        ],
        "modelReactionReagents": [
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00003_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00004_e0"
            },
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00010_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00011_e0"
            },
            {
                "coefficient": -1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00020_e0"
            },
            {
                "coefficient": 1,
                "modelcompound_ref": "~/modelcompounds/id/cpd00022_e0"
            }
        ],
        "modelcompartment_ref": "~/modelcompartments/id/e0",
        "name": "R00209",
        "probability": 0,
        "protons": 0,
        "reaction_ref": "~/template/reactions/id/rxn00154_e"
    }
    di = dict()

    rdr = Reaction.fromJSON(rd, di)
    print(rdr)
