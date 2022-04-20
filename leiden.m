function result = leiden(X, partition_type, params)
arguments
    X
    partition_type="RBC"
    params.resolution=1   %can be a vector of values to run
    params.n_iterations=-1
    params.max_comm_size=0
    params.rng_seed=42
    params.doRelabel=true
end

%TODO: support building adjacency from data in X

% tic
result = params;

disp('Performing Leidenalg (python) clustering...')

%add path to python script if needed...
Pypath = py.sys.path;
MLpath=string(path).split(';');
sctoolpath=MLpath(contains(MLpath,'scToolbox'));
if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
    insert(Pypath,int32(0),sctoolpath);
end

%X should contain the graph from umap... need to convert to a csr matrix
%can't pass sparse to python

N=size(X,1);
[sources,targets,weights]=find(X);
N = int32(N); %cast to int for python
sources = int32(sources-1); %convert to zero-based index as well
targets = int32(targets-1);

if length(params.resolution)==1

% def find_partition(N,sources,targets,
% weights=None, node_sizes=None, initial_membership=None, partition_type="RBC",
% resolution=1, n_iterations=-1, max_comm_size=0, rng_seed=None)

kwargs=pyargs('weights', weights, 'resolution', params.resolution,'n_iterations', ...
    int32(params.n_iterations), 'max_comm_size', int32(params.max_comm_size), ...
    'rng_seed', int32(params.rng_seed));

partition=py.leiden.find_partition(N,sources,targets,partition_type, kwargs);

groups = py.numpy.array(partition.membership);
clusterID=double(groups);
clusterID=clusterID+1; %python is 0-based
K=max(clusterID); 

%relabel to be in decreasing size
nT=arrayfun(@(x)nnz(clusterID==x),1:K);
if params.doRelabel
    [nT,ix]=sort(nT,'descend');
    IDold=clusterID;
    for i=1:K
        clusterID(IDold==ix(i))=i;
    end
end

else
end


% result.partition=partition;
result.K=K;
result.clusterID=clusterID;
result.mod=partition.modularity;
result.q=partition.q;
result.clusterCounts = nT;