import os, time
import numpy as np
import pandas as pd

from kasearch import AlignSequences, SearchDB
from PLAbDab.util import add_url_to_paired_data


minimum_size_of_results = 1_000_000


class SequenceSearch:
    """
    The main class for searching PLAbDab with KA-Search for sequences most similar to the query. 
    ...
    Methods
    -------
    sequence_search()
        Search PLAbDab with KA-Search for sequences most similar to the query.
    """

    def sequence_search(
        self, 
        seqs: dict, 
        keep_best_n: int = 10, 
        seq_identity_cutoff: float = 0.0,
        sort_by: str = 'average',
        regions=['whole'], 
        length_matched=[False],  
        allowed_species=['Human', 'Mouse'],
        url = False,
        n_jobs=5,
    ):
        """Search database with KA-Search for sequences most similar to the queries.
        
        Parameters
        ----------
        seqs : dict of str
            Heavy and light chain to search with. If only one chain, search unpaired data.
        keep_best_n : int
            Number of closest matches to return (default is 10).
        seq_identity_cutoff : float
            Lowest sequence identity returned (default is 0.0).
        sort_by : str
            Sort returned sequences by heavy, light or average (default is average).
        regions : list of str
            Region is search across (default is ['whole']).
        length_matched : list of bool
            If search only returns regions with the same length as query (default is [False]).  
        allowed_species : list of str
            Which species to search against (default is ['Human', 'Mouse']).
        url : bool
            If return results with additional column containing the url for the data (default is False)
        n_jobs : int
            Number of threads to use.
            
        """ 
        assert sort_by in ['heavy', 'light', 'average'], f"{sort_by} is invalid, use 'heavy', 'light' or 'average'"
        
        if len(seqs) == 1: 
            seq = seqs["H"] if "H" in seqs else seqs["L"]
            return search_unpaired(
                seq, 
                os.path.join(self.path_to_db, "kasearch_db"), 
                regions=regions, 
                length_matched=length_matched,  
                allowed_species=allowed_species, 
                keep_best_n=keep_best_n, 
                seq_identity_cutoff=seq_identity_cutoff,
                n_jobs=n_jobs,
            )
        
        elif len(seqs) == 2:
            out = search_paired(
                seqs["H"], 
                seqs["L"], 
                os.path.join(self.path_to_db, "kasearch_db"), 
                regions=regions, 
                length_matched=length_matched,  
                allowed_species=allowed_species, 
                keep_best_n=keep_best_n, 
                seq_identity_cutoff=seq_identity_cutoff,
                n_jobs=n_jobs,
                sort_by=sort_by,
            )
            if url:
                return add_url_to_paired_data(out)
            else:
                return out
        else:
            assert False, f"{seqs} needs to be a dict of either the heavy and light chain, or a single chain."


def search_unpaired(
    seq,
    db_path, 
    regions=['whole'], 
    length_matched=[False], 
    allowed_species=['Human', 'Mouse'], 
    keep_best_n=10, 
    seq_identity_cutoff=0.0,
    n_jobs=1,
):
    
    aligned_seqs = AlignSequences(allowed_species=allowed_species, n_jobs=n_jobs)(seq)
    
    chain_db = SearchDB(database_path=db_path, allowed_chain='Mix', regions=regions, length_matched=length_matched)
    chain_db.search(aligned_seqs[:1], keep_best_n=keep_best_n)
    output = chain_db.get_meta(n_query = 0, n_region = 0, n_sequences = 'all', n_jobs = n_jobs)

    return output[output.Identity >= seq_identity_cutoff].copy()


def search_paired(
    heavy_seq, 
    light_seq, 
    db_path, 
    sort_by = 'average', # heavy, light or average
    regions=['whole'], 
    length_matched=[False], 
    allowed_species=['Human', 'Mouse'], 
    keep_best_n=10, 
    seq_identity_cutoff=0.0,
    n_jobs=1,
):
    
    aligned_seqs = AlignSequences(allowed_species=allowed_species, n_jobs=n_jobs)([heavy_seq, light_seq])

    heavy_db = SearchDB(database_path=db_path, allowed_chain='Heavy', regions=regions, length_matched=length_matched)
    light_db = SearchDB(database_path=db_path, allowed_chain='Light', regions=regions, length_matched=length_matched)
    
    heavy_db.search(aligned_seqs[:1], keep_best_n=minimum_size_of_results)
    light_db.search(aligned_seqs[1:], keep_best_n=minimum_size_of_results)

    return find_best_pairs(heavy_db, light_db, keep_best_n=keep_best_n, seq_identity_cutoff=seq_identity_cutoff, sort_by=sort_by, n_jobs=n_jobs)


def find_best_pairs(
    heavy_db, 
    light_db, 
    keep_best_n=10, 
    sort_by='average',
    seq_identity_cutoff=0.0,
    n_jobs = 1
):

    heavy = np.concatenate([heavy_db.current_best_identities[0], heavy_db.current_best_ids[0,:,0,1:]], axis=1)
    light = np.concatenate([light_db.current_best_identities[0], light_db.current_best_ids[0,:,0,1:]], axis=1)
    heavy = pd.DataFrame(heavy, columns=['identity_heavy', 'idxs'])
    light = pd.DataFrame(light, columns=['identity_light', 'idxs'])

    paired = heavy.merge(right=light, left_on='idxs', right_on='idxs')
    paired.eval("identity_average = (identity_heavy + identity_light)/2", inplace=True)
    
    keep_best_n = min(len(paired) - 1, keep_best_n) 
    paired = paired.sort_values(by=f'identity_{sort_by}', ascending=False)[:keep_best_n].reset_index(drop=True)

    meta_idx = np.zeros((len(paired),2), dtype=np.int32)
    meta_idx[:,1] = paired.idxs.values.astype(np.int32)
    meta_data = heavy_db._extract_meta(meta_idx, n_jobs=n_jobs)

    paired_results = pd.concat([meta_data, paired], axis=1)
    paired_results.drop(columns=['Species','Chain', 'idxs'], inplace=True, errors='ignore')
    
    return paired_results[paired_results[f'identity_{sort_by}'] >= seq_identity_cutoff].copy()
    

