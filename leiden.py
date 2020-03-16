import leidenalg as la
import igraph as ig
import numpy as np
import time
from sklearn.neighbors import kneighbors_graph

# read in PCA coords from CSV
# coords = np.loadtxt('../matlab/analysis_Patient_3aggr/cleaned_pca_coords.csv',delimiter=',')
# k_NN = int(np.sqrt(coords.shape[0]))

# generate a KNN graph for Leiden
#  UMAP
# from umap.umap_ import nearest_neighbors
# from sklearn.utils import check_random_state
# random_state = check_random_state(0)

# start_time = time.time()
# knn_indices, knn_dists, forest = nearest_neighbors(X=coords, n_neighbors=k_NN, metric='euclidean',metric_kwds={},angular=False,random_state=random_state)
# print(f"umap nearest_neighbors took {time.time()-start_time}s")

def run(coords, k_NN, resolution):

    k_NN = int(k_NN) 

    # sklearn knn graph on UMAP embedding
    kw_knn_args = {'n_neighbors': k_NN, 'mode': 'distance', 'n_jobs': -1}

    start_time = time.time()
    adjacency = kneighbors_graph(coords, **kw_knn_args)
    adjacency.sort_indices()
    adjacency.eliminate_zeros()
    vcount = max(adjacency.shape)
    print(f"kneighbors_graph took {time.time()-start_time}s")

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    weights = np.transpose(weights)

    # now run leiden
    kw_la_args = {
        'n_iterations': 4,
        'seed': 42
    }
    kw_la_args['weights'] = weights
    # kw_la_args['partition_type'] = la.CPMVertexPartition
    # kw_la_args['resolution_parameter'] = resolution
    kw_la_args['partition_type'] = la.RBConfigurationVertexPartition
    kw_la_args['resolution_parameter'] = resolution

    edgelist = list(zip(sources.tolist(), targets.tolist()))
    G = ig.Graph(vcount, edgelist)
    start_time = time.time()
    partition = la.find_partition(G, **kw_la_args)
    labels = np.array(partition.membership)
    nclust = np.unique(labels).max()+1
    print(f"leidenalg took {time.time()-start_time}s")
    print(f"number of clusters = {nclust}")

    return labels
    # write the labels to CSV


# plot
# fig = plt.figure(figsize=(4, 4), dpi=144)
# ax = fig.add_subplot(111)
# palette = sns.color_palette('deep', np.unique(labels).max() + 1)
# colors = [palette[x] for x in labels]
# plt.scatter(coords[:, 0], coords[:, 1], c=colors)
# plt.title(f'{nclust} clusters found by Leiden', fontsize=8)
# frame = plt.gca()
# frame.axes.get_xaxis().set_visible(False)
# frame.axes.get_yaxis().set_visible(False)
# plt.tight_layout()
# # plt.savefig(figfileroot + 'umap_modularity.png')
# plt.show()