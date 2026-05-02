import itertools
from Bio.SeqUtils import molecular_weight
import gzip
import csv
import sys
import pandas as pa
from time import time

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
pio.templates.default = "plotly_white"

from pathlib import Path
project_root = str(Path(__file__).resolve()).split('src')[0]
sys.path.append(project_root)
from src.util.modelComponents import *
from src.util.bcolors import bcolors

# Atomic mass of Fe is 55.845 g/mol
# Atomic mass of Zn is 65.38 g/mol

avogadro = 6.02214076e+23

class ProteinWeightGenerator:
    def __init__(self, spc_name:str="Sorghum", project='QPSI', rnaseq_path:str=None, model_genes:set={}, project_cols:list=[], cap_percent:int=-1):

        self.spc_name               = spc_name
        self.model_genes            = model_genes
        self.rnaseq_path            = rnaseq_path
        self.project                = project
        self.project_cols           = project_cols if project_cols else ['tissue', 'treatment', 'time_stamp']
        self.cap_percent            = cap_percent

        ## Plastid proteome data 
        self.chloro_atha_file = os.path.join(project_root, "data", "plastid_proteome", "Full_AT_CHLORO_2019.tsv")
        self.PPDB_file = os.path.join(project_root, "data", "plastid_proteome", "PPDB_Plastid_Proteome.tsv")

        ## Amino Acids propeties including "molecular_mass" and "num_atoms"
        self.AA_properties_file = os.path.join(project_root, "data", "ProteinSeq", "AA_properties.csv")

        ## Protein seqeuences: A. Thaliana and other species
        ## And orthologs files     
        self.protein_weight_file    = os.path.join(project_root, "data", "ProteinSeq", f"{self.spc_name}_protein_weights.csv")
        ortho_version = 'apr10'
        if spc_name.lower() == 'sorghum':
            self.spc_fasta          = os.path.join(project_root, "data", "ProteinSeq", "Sbicolor_454_v3.1.1.protein.fa")
            self.orthologs_file     = os.path.join(project_root, "data", "orthologs", f"orthologs_{ortho_version}", f"Sbicolor_v3.1.1_Athaliana_Araport11_Functional_Homologs_{ortho_version}.tsv")

        elif spc_name.lower() == 'poplar':
            self.spc_fasta          = os.path.join(project_root, "data", "ProteinSeq", "Ptrichocarpa_533_v4.1.protein.fa")
            self.orthologs_file     = os.path.join(project_root, "data", "orthologs", f"orthologs_{ortho_version}", f"Ptrichocarpa_v4.1_Athaliana_Araport11_Functional_Homologs_{ortho_version}.tsv")
            
        else:
            self.spc_fasta          = os.path.join(project_root, "data", "ProteinSeq", "Athaliana_447_Araport11.protein.fa")
            self.orthologs_file     = None

        ## for data processing
        self.ortho_mapping_dict     = dict()
        self.protein_seqs_dict      = dict()
        self.AA_atoms_dict          = dict()
        self.protein_weight_dict    = dict()
        self.protein_moles_dict     = dict()
        self.atha_chloro_proteins   = set()
        self.plastid_weight_sums    = None

    def readAthatChloroProteins(self):
        print(bcolors.PROG+f"Reading A. thaliana chloroplast proteins ..."+bcolors.ENDC)

        with open(self.chloro_atha_file) as fh:
            reader = csv.reader(fh, delimiter='\t')
            # This skips the first row of the CSV file.
            next(reader)

            for row in reader:
                if not row:
                    continue

                if(row[0]==''):
                    continue

                if('?' in row[3]):
                    continue

                if('ENV' not in row[3]) and ('THY' not in row[3]) and ('STR' not in row[3]):
                    continue

                self.atha_chloro_proteins.add(row[1])

        with open(self.PPDB_file) as fh:
            reader = csv.reader(fh, delimiter='\t')
            # This skips the first row of the CSV file.
            next(reader)

            for row in reader:
                if not row:
                    continue

                if(row[0]==''):
                    continue

                self.atha_chloro_proteins.add(row[0])

    # load RNASeq data (transcriptome) to a pandas dataframe
    def readRNASeq(self):
        RNASeq_file_path = self.rnaseq_path

        rnaseq_df = pa.read_csv(RNASeq_file_path)
        rnaseq_df = rnaseq_df[["Gene_ID", "value"]+self.project_cols]

        if ('tmm' in self.rnaseq_path) and (self.project == 'QPSI'): 
            if (self.cap_percent > 0):
                print(bcolors.PROG+"ProteinWeight Generator: Value cap at ", self.cap_percent, bcolors.ENDC)
                rnaseq_df.loc[(rnaseq_df["value"] > self.cap_percent), "value"] = self.cap_percent

        return rnaseq_df

    def readOrthologs(self, files=[], verbose=False):
        if not self.atha_chloro_proteins:
            self.readAthatChloroProteins()
        
        if verbose: print(bcolors.PROG+"There are ", len(self.atha_chloro_proteins), " A. Thaliana platid proteins"+bcolors.ENDC)

        if(len(files)==0):
            files.append(self.orthologs_file)

        for file in files:
            print(bcolors.PROG+f"Generating orthologs dict from {file} ..."+bcolors.ENDC)

            with open(file, mode='r') as inFile:
                reader = csv.reader(inFile, delimiter='\t')
                for row in reader:
                    AT_genes = {gene.rsplit('.', 1)[0] for gene in row[0].split(', ')}
                    for at_gene in AT_genes:
                        # keep proteins found in the chloroplast only
                        if at_gene in self.atha_chloro_proteins:
                            spc_genes = {gene.rsplit('.', 1)[0] for gene in row[1].split(', ')}
                            try:
                                self.ortho_mapping_dict[at_gene].update(spc_genes)
                            except KeyError:
                                self.ortho_mapping_dict[at_gene] = spc_genes

        if not self.ortho_mapping_dict:
            print(bcolors.WARNING+"No ortholog mapping found!"+bcolors.ENDC)
            sys.exit(1)

    # Read protein (RNA) sequences from provided fasta files
    def readProteinSeq(self, chloro_only=False, model_genes_only=False, prim_seq=False, file_configs=[]):
        # All species orthologs found in the plastid
        chloro_proteins = set()
        if (not self.ortho_mapping_dict) and chloro_only:
            if 'atha' not in self.spc_name.lower():
                self.ortho_mapping_dict = self.readOrthologs_new()
                chloro_proteins = set().union(*map(set, [v for v in self.ortho_mapping_dict.values()]))

        print(bcolors.PROG+f"Reading protein sequences ..."+bcolors.ENDC)

        if(len(file_configs)==0):
            file_configs.append({'file':self.spc_fasta,'gene_regex':r'^(\S+)'})

        for config in file_configs:
            file = config['file']
            gene_id_pattern = re.compile(config['gene_regex'])

            # open fasta file and alternate header and sequence
            if self.is_gzipped(file):
                fasta_handle = gzip.open(file, 'rt')
            else:
                fasta_handle = open(file,'r')
            fasta_iterator = (x[1] for x in itertools.groupby(fasta_handle, lambda line: line[0] == ">"))
            for header in fasta_iterator:
                # drop the ">"
                header = header.__next__()[1:].strip()
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in fasta_iterator.__next__())
                # extract the ID field from header
                gene_id = ProteinWeightGenerator.getGeneID(header, gene_id_pattern)

                # check if gene is in chloroplast if chloro_only==True
                if chloro_only and (gene_id not in chloro_proteins):
                    continue
                # process model gene only if model_genes_only == True
                if model_genes_only and self.model_genes and (gene_id not in self.model_genes):
                    continue

                # Associate multiple sequences with a single ID by appending to a list
                self.protein_seqs_dict.setdefault(gene_id, []).append(seq.rstrip('*'))

                # Compute protein molecular weight
                self.protein_weight_dict.setdefault(gene_id, []).append(molecular_weight(seq.rstrip('*'),"protein"))

        print(bcolors.PROG+f"  -> There are {len(self.protein_seqs_dict)} sequences ..."+bcolors.ENDC)

    # Process the header file to extract the gene ID and the transcripy ID
    # return gene_ID and transcript ID
    @staticmethod
    def getGeneID(line, id_pattern=None):
        if id_pattern:
            match = id_pattern.search(line)
            if match:
                return match.group(1)

        # Fallback: use the first word of the header text if no pattern or match fails
        # Safely drops trailing description metadata
        return line.split()[0]

    # Read amino acids properites form file "AA_properties_file"
    # returns a dictionary with "abreviation_s"->"num_atoms"
    def read_AA_properties(self):
        sep = ','

        with open(self.AA_properties_file, mode='r') as inFile:
            reader = csv.reader(inFile, delimiter=sep)
            for row in reader:
                # skip header
                if row[0] == 'AA_name':
                    continue

                # Build AA atoms dictionary
                self.AA_atoms_dict[row[2]] = int(row[4]) if row[4] else 0

        if not self.AA_atoms_dict:
            sys.exit(1)

    # Compute the mass based on the sequence
    # using the "AA_atoms_dict"
    def compute_protein_atoms(self, seq):
        num_atoms = 0

        for AA in self.AA_atoms_dict.keys():
            num_atoms += seq.count(AA) * self.AA_atoms_dict[AA]

        return num_atoms/avogadro

    # Compute the protein weight and save it in "self.protein_weight_dict"
    # compute the total plastid protein mass and save it in "self.plastid_weight_sums"
    def computeWeightsAndMass(self, to_file=False):
        totalMass = 0
        rnaseq_df = self.readRNASeq()
        rnaseq_genes = set(rnaseq_df["Gene_ID"].unique())

        # Read Amino Acids molecular mass
        if not self.AA_atoms_dict:
            self.read_AA_properties()

        # Read chloroplast DB
        if not self.atha_chloro_proteins:
            self.readAthatChloroProteins()

        # Read Athaliana to Species otholog mapping
        if not self.ortho_mapping_dict:
            self.ortho_mapping_dict = self.readOrthologs_new(mode=2)

        # Read PPDB plastid proteome
        PPDB_df = pa.read_csv(self.PPDB_file, sep='\t')
        PPDB_df['Accession'] = PPDB_df['Accession'].apply(lambda x: x.split('.')[0] )
        PPDB_genes = set(PPDB_df['Accession'].unique())

        # List of all Arabidopsis genes that are in the both lists:
        #    (atha_chloro_proteins and PPDB_genes)
        atha_plastid_genes = PPDB_genes.union(self.atha_chloro_proteins)
        # Find plastidial species genes using athaliana plastidial genes orthologs
        spc_plastid_genes  = {gene for at_gene in atha_plastid_genes
                                    if at_gene in self.ortho_mapping_dict
                                    for gene in self.ortho_mapping_dict[at_gene]}

        ## printing information about the data
        in_plastid_not_inRNA = spc_plastid_genes-rnaseq_genes
        print('Atha plastid genes: ', len(atha_plastid_genes))
        print("Num.", self.spc_name, " plastidial genes: ", len(spc_plastid_genes), \
            " num. of genes in the RNASeq file: ", len(rnaseq_genes), \
            " -> genes in the plastid but not present in the RNASeq file: ",\
            len(in_plastid_not_inRNA), '\n', in_plastid_not_inRNA)
        
        # Read protein sequences to compute the protein weight
        if not self.protein_seqs_dict:
            self.readProteinSeq()

        # Compute the number of atoms in each protein
        if not self.protein_moles_dict:
            self.protein_moles_dict  = {gene_id:self.compute_protein_atoms(seq)
                                            for gene_id, seq in self.protein_seqs_dict.items()}

        #------- Compute protein moleculare weight
        self.protein_weight_dict = {gene_id:molecular_weight(seq,"protein")
                                        for gene_id, seq in self.protein_seqs_dict.items()}

        # write results to file
        if to_file:
            with open(self.protein_weight_file, 'w') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(['Gene_ID', "weight", "moles", "g_mass"])
                sum = 0
                for key, value in self.protein_weight_dict.items():
                    moles = self.protein_moles_dict[key]
                    sum += moles*value
                    writer.writerow([key, value, moles, moles*value])

        #------- Compute total plastid weight
        # Compute protein weight using the "protein_moles_dict"
        rnaseq_df = rnaseq_df[rnaseq_df['Gene_ID'].isin(spc_plastid_genes)]
        
        if self.project == 'QPSI':
            rnaseq_df['time_stamp'] = rnaseq_df['time_stamp'].transform(lambda ts:ts if ts in ['14d', '21d'] else '0'+ts)
        rnaseq_df['weight'] = rnaseq_df['Gene_ID'].map(self.protein_moles_dict)
        rnaseq_df['weight'] = rnaseq_df['weight']#*avogadro
        rnaseq_df["TotalPlastidMass"] = rnaseq_df['weight']*rnaseq_df['value']

        # compute total plastid protein mass for each col "project_cols"
        weight_sums =  rnaseq_df.groupby(self.project_cols)\
        .agg({'TotalPlastidMass':'sum'}).reset_index()
        
        if to_file: weight_sums.to_csv(f'plastidProteinMass{self.spc_name}.tsv', sep='\t', index=False)
        
        weight_sums.set_index(self.project_cols, inplace=True)
        self.plastid_weight_sums = weight_sums
        return weight_sums

    # Plot the protein mass across different treatments, time points, tissues, etc.
    def plot_plastid_weight(weight_sums):
        # This is for QPSI data
        dss = ['02d', '04d', '07d', '14d', '21d']
        weight_sums = weight_sums[weight_sums['time_stamp'].isin(dss)]
        print(weight_sums)
        fig = px.line(weight_sums, x='time_stamp', y="TotalPlastidMass", color='treatment', facet_row='tissue', category_orders={"time_stamp": dss}, width=600, title=f"{self.spc_name}")

        fig.show()

    @staticmethod
    def is_gzipped(filepath):
        # Peek at the first two bytes in binary mode
        with open(filepath, 'rb') as test_f:
            return test_f.read(2) == b'\x1f\x8b'
    
    ## Additional processing for gene localization and function
    @staticmethod
    def process_localization(ftr_loc_dict, loc_dict, role='', reactions=[], publications=[]):
        for cpt, data in loc_dict.items():
            # Example data item: 
            # "d": {
            #     "Athaliana_TAIR10||AT2G14750": [
            #         "PPDB"
            #     ],
            #     "Athaliana_TAIR10||AT4G39940": [
            #         "PPDB"
            #     ]
            # },
            for gene, loc in data.items():
                if '||' in gene:
                    gene = gene.split('||')[1]
                    if gene in ftr_loc_dict:
                        ftr_loc_dict[gene].append({'localization': cpt,
                                              'source' : loc,
                                              'role': role,
                                              'reactions': reactions,
                                              'publications': publications})
                    else:
                        ftr_loc_dict[gene] = [{'localization': cpt,
                                              'source' : loc,
                                              'role': role,
                                              'reactions': reactions,
                                              'publications': publications}]
        return ftr_loc_dict


if __name__ == "__main__":
    spc = "Poplar"
    RNAseq_path = os.path.join(project_root, "data", "RNASeq_data", "tmm", f"{spc}_raw_genes_tmm_mean_std.csv")

    pwgObj = ProteinWeightGenerator(spc_name=spc, rnaseq_path=RNAseq_path)
    pwgObj.computeWeightsAndMass(to_file=True)
    
 
