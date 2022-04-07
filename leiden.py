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


def get_adjacency(coords, k_NN=None):
    if k_NN is None:
        k_NN = np.sqrt(coords.shape[0])

    k_NN = int(k_NN)

    # sklearn knn distance graph... would be better to do SNN?
    kw_knn_args = {'n_neighbors': k_NN, 'mode': 'distance', 'n_jobs': -1}

    start_time = time.time()
    adjacency = kneighbors_graph(coords, **kw_knn_args)
    adjacency.sort_indices()
    adjacency.eliminate_zeros()
    print(f"kneighbors_graph took {time.time()-start_time}s")

    return adjacency

# def run(adjacency, resolution, k_NN=None, rng_seed=42, n_iterations=4):
#     N = adjacency.shape[0]
#     # build igraph object
#     sources, targets = adjacency.nonzero()
#     weights = adjacency[sources, targets]
#     if isinstance(weights, np.matrix):
#         weights = weights.A1


def find_partition(
    N,
    sources,
    targets,
    partition_type,
    weights=None,
    node_sizes=None,
    initial_membership=None,
    resolution=1,
    n_iterations=-1,
    max_comm_size=0,
    rng_seed=None
):

    # graph from adjacency matrix (passed in as sparse representation: N, sources, targets, weights)
    g = ig.Graph(directed=True)
    g.add_vertices(N)  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))

    # handle the possible partitions and their allowed args
    partition_kwargs = {}
    if partition_type is None:
        partition_type = la.RBConfigurationVertexPartition

    elif partition_type == "Modularity":
        partition_type = la.ModularityVertexPartition
        resolution = None
        node_sizes = None

    elif partition_type == "RBC":
        partition_type = la.RBConfigurationVertexPartition
        node_sizes = None

    elif partition_type == "RBERVertexPartition":
        partition_type = la.RBERVertexPartition

    elif partition_type == "CPM":
        partition_type = la.CPMVertexPartition

    elif partition_type == "Significance":
        partition_type = la.SignificanceVertexPartition
        resolution = None
        weights = None

    elif partition_type == "Surprise":
        partition_type = la.SignificanceVertexPartition
        resolution = None

    else:
        print(
            f"bad partition_type: {partition_type}. Using RBConfigurationVertexPartition")
        partition_type = la.RBConfigurationVertexPartition

    if weights is not None:
        g.es['weight'] = weights

    if resolution is not None:
        partition_kwargs['resolution_parameter'] = resolution

    if node_sizes is not None:
        partition_kwargs['node_sizes'] = node_sizes

    start_time = time.time()

    # find a partition
    ## leidenalg.find_partition(graph, partition_type, initial_membership=None, weights=None, n_iterations=2, max_comm_size=0, seed=None, **kwargs)
    partition = la.find_partition(g, partition_type, initial_membership=initial_membership, weights=weights,
                                  n_iterations=n_iterations, max_comm_size=max_comm_size, seed=rng_seed, **partition_kwargs)

    groups = np.array(partition.membership)
    print(f"leidenalg took {time.time()-start_time}s")
    print(f"number of clusters = {len(np.unique(groups))}")

    return partition
    # write the labels to CSV


def resolution_profile(
    N,
    sources,
    targets,
    partition_type,
    resolution_range,
    weights=None,
    min_diff_bisect_value=1,
    min_diff_resolution=0.001,
    linear_bisection=False,
    number_iterations=1,
    rng_seed=None
):

    # graph from adjacency matrix (passed in as sparse representation: N, sources, targets, weights)
    g = ig.Graph(directed=True)
    g.add_vertices(N)  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))

    # handle the possible partitions and their allowed args
    partition_kwargs = {}
    if partition_type is None:
        partition_type = la.RBConfigurationVertexPartition

    elif partition_type == "RBC":
        partition_type = la.RBConfigurationVertexPartition
        node_sizes = None

    elif partition_type == "RBERVertexPartition":
        partition_type = la.RBERVertexPartition

    elif partition_type == "CPM":
        partition_type = la.CPMVertexPartition

    else:
        print(
            f"bad partition_type: {partition_type}. Using RBConfigurationVertexPartition")
        partition_type = la.RBConfigurationVertexPartition

    if weights is not None:
        g.es['weight'] = weights

    if node_sizes is not None:
        partition_kwargs['node_sizes'] = node_sizes

    start_time = time.time()

    optimiser = la.Optimiser()

    if rng_seed is not None:
        optimiser.set_rng_seed(rng_seed)

    # resolution_profile(graph, partition_type, resolution_range, weights=None, bisect_func=<function Optimiser.<lambda>>,
    #                   min_diff_bisect_value=1, min_diff_resolution=0.001, linear_bisection=False, number_iterations=1, **kwargs)
    profile = optimiser.resolution_profile(g, partition_type, resolution_range, weights=weights, min_diff_bisect_value=min_diff_bisect_value, 
        min_diff_resolution=min_diff_resolution, linear_bisection=linear_bisection, number_iterations=number_iterations, **partition_kwargs)

    print(f"leidenalg took {time.time()-start_time}s")
    # print(f"number of clusters = {len(np.unique(groups))}")

    return profile
