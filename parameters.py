from pathlib import Path
import os
import sys

project_root = str(Path(__file__).resolve().parent)
sys.path.append(project_root)
runtime_root = os.getcwd()

# Global Species Synonym Mappings
parameters_all_spc = {
    "TSU": {'name': 'TSU', 'synonyms': ['TSU', 'Athaliana']},
    "C24": {'name': 'C24', 'synonyms': ['C24', 'Athaliana']},
    "Athaliana": {'name': "Arabidopsis", 'synonyms': ["Atha", "athaliana", "Arabidopsis", 'sandbox', 'Athaliana']},
    "Poplar": {'name': "Poplar", 'synonyms': ["Poplar", "Pptrich", "Ptrichocarpa", "ptrich", "ptr"]},
    "Sorghum": {'name': "Sorghum", 'synonyms': ["Sorghum", "Sbicolor", "Sbi", "sbi"]}
}

class Parameters:
    def __init__(self, verbose=True):
        self.project         = 'qpsi-plastidial'
        self.value_column    = 'value'
        self.project_species = ["Sorghum", "Poplar"]
        self.group_columns   = ['condition']
        self.rnaSeq_id_col   = 'Gene_ID'
        self.error           = False
        self.control_id      = 'Control'
        self.control_name    = 'Control'  # Added to align with downstream scripts
        self.trmt_colmn      = 'condition'
        self.reaction_score_method = 'sum'
        self.ignore_organellar_roles = os.path.join(project_root, "data", "organellar-encoded_subunits_to_ignore.txt")

        # Baseline Directory Structuring
        self.json_files_folder = os.path.join(project_root, "projects", self.project, "inputs")
        self.RNASeq_folder     = os.path.join(project_root, "projects", self.project, "rnaseq-data")
        self.results_folder    = os.environ.get("RESULTS_FOLDER", os.path.join(project_root, "projects", self.project, "integration_results"))
        self.reaction_file_suffix = '_reaction_scores.tsv'

        # Initial Directory Existence Testing
        for folder_name, folder_path in [("Models Input", self.json_files_folder), 
                                         ("RNASeq Data", self.RNASeq_folder), 
                                         ("Results Output", self.results_folder)]:
            if not os.path.exists(folder_path):
                print(f"[!] WARNING: {folder_name} folder not found. Expected location: {folder_path}")
                if folder_name == "Results Output":
                    os.makedirs(folder_path, exist_ok=True)
                    print(f"    -> Automatically initialized folder path.")
                else:
                    self.error = True

        # Automated File Discovery & Integrity Mapping
        self.model_paths      = {}
        self.rnaseq_paths     = {}
        self.rxn_scores_paths = {}  # Added to house optional results files

        if not self.error:
            self._discover_json_models(verbose)
            self._discover_rnaseq_files(verbose)
            self._discover_existing_results(verbose)
            self._validate_discovered_assets()

    def _discover_json_models(self, verbose):
        """Scans the inputs folder to map active reconstructions via target synonyms."""
        if verbose: print("-" * 10 + " Updating JSON file paths " + "-" * 10)
        try:
            file_names = [os.path.join(self.json_files_folder, x) for x in os.listdir(self.json_files_folder)]
        except Exception as e:
            print(f" [FAIL] Scanning model repository failed: {e}")
            self.error = True
            return

        for spc in self.project_species:
            synonyms = parameters_all_spc.get(spc, {}).get('synonyms', [spc])
            for file_path in file_names:
                if "media" in os.path.basename(file_path).lower():
                    continue
                if any(y.lower() in os.path.basename(file_path).lower() for y in synonyms) and file_path.endswith('.json'):
                    self.model_paths[spc] = file_path
                    if verbose: print(f" -> Found {spc} Model: {os.path.basename(file_path)}")
                    break

    def _discover_rnaseq_files(self, verbose):
        """Scans specified resolution directories to find relevant expression tables."""
        if verbose: print("-" * 10 + " Updating RNASeq file paths " + "-" * 10)
        target_dir = self.RNASeq_folder
        if not os.path.exists(target_dir):
            print(f" [FAIL] Target RNASeq data subdirectory missing: {target_dir}")
            self.error = True
            return

        try:
            file_names = os.listdir(target_dir)
        except Exception as e:
            print(f" [FAIL] Reading target data path failed: {e}")
            self.error = True
            return

        for file_name in file_names:
            for spc in self.project_species:
                synonyms = parameters_all_spc.get(spc, {}).get('synonyms', [spc])
                if any(y.lower() in file_name.lower() for y in synonyms):
                    self.rnaseq_paths[spc] = os.path.join(target_dir, file_name)
                    if verbose: print(f" -> Found {spc} RNASeq Profile: {file_name}")
                    break

    def _discover_existing_results(self, verbose):
        """Looks strictly for files matching the standardized layout: {Species}{Suffix}"""
        if not os.path.exists(self.results_folder):
            return

        try:
            result_files = os.listdir(self.results_folder)
        except Exception as e:
            print(f" [FAIL] Reading results directory failed: {e}")
            return

        found_any_path = False
        for spc in self.project_species:
            synonyms = parameters_all_spc.get(spc, {}).get('synonyms', [spc])
            found_path = None
            
            # Formulate the explicit target filename syntax using your suffix property
            # matching layouts like "Sorghum_reaction_score.tsv" or "Sbicolor_reaction_score.tsv"
            for file_name in result_files:
                if file_name.endswith(self.reaction_file_suffix):
                    # Extract the prefix to see if it belongs to this species mapping group
                    base_prefix = file_name[:-len(self.reaction_file_suffix)]
                    if any(y.lower() == base_prefix.lower() for y in synonyms):
                        found_path = os.path.join(self.results_folder, file_name)
                        break
            
            if found_path:
                self.rxn_scores_paths[spc] = found_path
                found_any_path = True
            else:
                self.rxn_scores_paths[spc] = None
                # if verbose: print(f" -> Optional: No files generated yet matching prefix for {spc} (Targeting: *{self.reaction_file_suffix})")
        if(found_any_path):
            if verbose: print("-" * 10 + " Checked for Integrated Results files " + "-" * 10)
            for spc in self.rxn_scores_paths:
                if verbose: print(f" -> Found {spc} Results File: {os.path.basename(self.rxn_scores_paths[spc])}")


    def _validate_discovered_assets(self):
        """Verifies all explicitly targeted workspace constraints were mapped cleanly."""
        for spc in self.project_species:
            if spc not in self.model_paths:
                print(f" [ERR] Missing JSON construction asset! No valid model found for target species: {spc}")
                self.error = True
            if spc not in self.rnaseq_paths:
                print(f" [ERR] Missing Core Transcript Balance Metrics! No RNASeq profile found matching target species: {spc}")
                self.error = True

        if self.error:
            print("\n [CRITICAL FAILURE] Project workspace integrity checks failed. Review error reports.")


if __name__ == "__main__":
    param = Parameters(verbose=True)
    if not param.error:
        print("\n=== All Systematic Checks Passed! ===")
        print(f"Results Trackings: {param.rxn_scores_paths}")