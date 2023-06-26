import os 
import pandas as pd
from ImmuneBuilder import refine
import PLAbDab.structure_search_utils as utils
from PLAbDab.util import add_url_to_paired_data

class StructureSearch:
    
    
    def __init__(self, path_to_db):
        
        self.path_to_db = path_to_db
        
        self.paired_sequences = pd.read_csv(os.path.join(self.path_to_db, "paired_sequences.csv.gz"))
        
        
    
    def structure_search(self, seqs, rmsd_cutoff = 10.0, url = False, save_query = False, filename = "temp_structure.pdb"):
        assert ("H" in seqs) and ("L" in seqs), f"{seqs} needs to be a dict of a heavy and light chain."

        # Get path to models
        model_db_path = os.path.join(self.path_to_db, 'models')

        # Model antibody
        filename = utils.fast_unrefined_antibody_model(seqs, filename=filename)
        antibody = utils.parse_antibody(filename)

        # If the model of the query should be saved, refine it. Otherwise, delete it
        if save_query:
            refine.refine(filename, filename)
        else:
            os.remove(filename)

        # Find CDR_lengths
        CDR_length = "_".join(utils.get_CDR_lengths(antibody))

        # Find antibodies with same CDR lengths
        cdr_length_cluster = self.paired_sequences[self.paired_sequences.cdr_lengths == CDR_length].copy()

        # If there are none we break
        if len(cdr_length_cluster) == 0:
            print("No structures with the same CDR lengths")
            return None

        # Take unique models (don't calculate rmsd twice on same structure)
        unique_cluster = cdr_length_cluster.drop_duplicates(["model"])
        
        # Iterate through entries with same length CDRs
        rmsds = {}
        for model_name in unique_cluster.model:
            if model_name in ["None", "FAILED"]:
                rmsds[model_name] = float("Nan")
                continue
            
            path = os.path.join(model_db_path, model_name + ".pdb")

            # If model does not exist skip it
            if os.path.exists(path):
                db_antibody = utils.parse_antibody(path)
                rmsds[model_name] = utils.rmsd(antibody, db_antibody)
            else:
                rmsds[model_name] = float("Nan")
        
        cdr_length_cluster["rmsd"] = [rmsds[model_name] for model_name in cdr_length_cluster.model]

        output = cdr_length_cluster.sort_values(by="rmsd").reset_index(drop=True)
        output = output[output.rmsd <= rmsd_cutoff]
        
        if url:
            output = add_url_to_paired_data(output)
        
        return output