import os
import json
import pandas as pd
import torch

from PLAbDab.sequence_search import SequenceSearch
from PLAbDab.structure_search import StructureSearch
from PLAbDab.util import add_url_to_paired_data

class PLAbDab(SequenceSearch, StructureSearch):
    """
    The main class for searching and retrieving data from the Patent and Literature Antibody Database (PLAbDab). 
    ...
    Attributes
    ----------
    path_to_db : str
        Path to local PLAbDab database (up to date version can be downloaded from https://opig.stats.ox.ac.uk/webapps/plabdab/).
    n_jobs : int
        Number of threads to use (default 5).
    Methods
    -------
    sequence_search()
        Search PLAbDab with KA-Search for sequences most similar to the query.
    structure_search()
        Search PLAbDab with ABB2 for structures most similar to the query.
    get_db()
        Retrieve data from PLAbDab.
    column_search()
        Search columns in PLAbDab for specific terms.
    """
    
    def __init__(
        self, 
        path_to_db, 
        n_jobs = 5, 
        **kwargs
    ):
        super().__init__(path_to_db)
        self.config = kwargs
        self.path_to_db = path_to_db
        self.__check_config()
        
        
        if not os.path.exists(os.path.join(self.path_to_db, "paired_sequences.csv.gz")):
            raise ValueError(f"Provided path_to_db ({self.path_to_db}) does not contain a paired_sequences.csv.gz file.")
            
        torch.set_num_threads(n_jobs)
        self.paired_sequences = pd.read_csv(os.path.join(self.path_to_db, "paired_sequences.csv.gz"))
        self.unpaired_sequences = pd.read_csv(os.path.join(self.path_to_db, "unpaired_sequences.csv.gz"))

    def __check_config(self):
        config_file = os.path.join(self.path_to_db, "config.json")

        if os.path.exists(config_file):
            with open(config_file, 'r') as fp:
                config = json.load(fp)

            for key in config:
                if key in self.config:
                    assert config[key] == self.config[key], f"Database was initiallised with a different configuration:\n{self.config}"
                else:
                    self.config[key] = config[key]

    def get_db(self, paired=True, url=True):

        if paired:
            if not url:
                return self.paired_sequences
            else:
                return add_url_to_paired_data(self.paired_sequences)
        else:
            return self.unpaired_sequences

    def column_search(self, term = "antibod", paired=True, column=None, url=True, case = False):
        if paired:
            df = self.paired_sequences
            if column is None:
                column = "reference_title"
        else:
            df = self.unpaired_sequences
            if column is None:
                column = "GBReference_title"

        output = df[df[column].str.contains(term, case=case, na=False)].reset_index(drop=True)
        if url and paired:
            output = add_url_to_paired_data(output)

        return output