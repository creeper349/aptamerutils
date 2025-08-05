import numpy as np
from tqdm import tqdm
from typing import Literal
from multiprocessing import Pool, cpu_count
from .sequence import Sequence
from itertools import combinations
from functools import partial
from ._posmap import PosMap

class DistMatrix(np.ndarray):
    def __new__(cls, input_array, labels=None):
        obj = np.asarray(input_array).view(cls)
        obj.labels = labels
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.labels = getattr(obj, 'labels', None)
        
    @classmethod
    def fromSequencesParallel(cls, seqs: list[Sequence], dist_func:Literal["ed", "ncd", "kmer"] = "kmer", verbose = True, **kwargs):
        n = len(seqs)
        compute = partial(compute_pair, seqs = seqs, dist_func = dist_func, **kwargs)
        
        matrix = np.zeros((n, n))
        pairs = list(combinations(range(n), 2))
        with Pool(processes = min(cpu_count(), 8)) as pool:
            results = list(tqdm(pool.imap(compute, pairs), total=len(pairs), disable=not verbose))
        
        for i, j, val in results:
            matrix[i, j] = matrix[j, i] = val
        return cls(matrix)
    
    @classmethod
    def fromSequences(cls, seqs: list[Sequence], dist_func:Literal["ed", "ncd", "kmer"] = "kmer", verbose = True, **kwargs):
        n = len(seqs)
        
        matrix = np.zeros((n, n))
        pairs = list(combinations(range(n), 2))
        
        for i, j in tqdm(pairs, total = len(pairs), disable= not verbose):
            matrix[i, j] = matrix[j, i] = compute_pair((i, j), seqs, dist_func, **kwargs)
        return cls(matrix)
    
    def createPosMap(self, mapping_metrics:Literal["umap", "mds"] = "umap", **kwargs):
        if mapping_metrics == "umap":
            matrix = PosMap.fromDistMapUMAP(self, **kwargs)
        else:
            matrix = PosMap.fromDistMapMDS(self, **kwargs)
        return matrix
    
def compute_pair(pair, seqs: list[Sequence], dist_func, **kwargs):
    i, j = pair
    return seqs[i].distWith(seqs[j], dist_func, **kwargs)
    
if __name__ == "__main__":
    pass