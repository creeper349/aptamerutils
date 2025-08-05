import numpy as np
from .sequence import Sequence
from typing import Literal
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from umap import UMAP

class PosMap(np.ndarray):
    def __new__(cls, input_array, labels=None):
        obj = np.asarray(input_array).view(cls)
        obj.labels = labels
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.labels = getattr(obj, 'labels', None)
        
    @classmethod
    def fromDistMapMDS(cls, dist, n_components = 20, verbose = True, **kwargs):
        mds = MDS(n_components, verbose = verbose, **kwargs)
        scaler = StandardScaler()
        return cls(scaler.fit_transform(mds.fit_transform(dist)))
    
    @classmethod
    def fromDistMapUMAP(cls, dist, n_neighbors = 50, n_components = 20, verbose = True, **kwargs):
        umap = UMAP(n_neighbors, n_components, metric = "precomputed", verbose = verbose, **kwargs)
        scaler = StandardScaler()
        return cls(scaler.fit_transform(umap.fit_transform(dist)))
    
    def getCluster(self, metrics:Literal["dbscan", "kmeans"] = "dbscan", **kwargs):
        if metrics == "dbscan":
            eps = kwargs.pop("eps", 0.6)
            dbscan = DBSCAN(eps, **kwargs)
            return dbscan.fit_predict(self)
        else:
            n_clusters = kwargs.pop("n_clusters", 8)
            kmeans = KMeans(n_clusters, **kwargs)
            return kmeans.fit_predict(self)