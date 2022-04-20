function result = leiden_profile(X, partition_type, resolution_range, params)
arguments
    X
    partition_type="RBC"
    resolution_range=[0,1]   %can be a vector of values to run
    params.min_diff_bisect_value=1,
    params.min_diff_resolution=0.001,
    params.linear_bisection=false,
    params.n_iterations=1,
    params.rng_seed=42
    params.doRelabel=true
end

%TODO: support building adjacency from data in X

% tic
result = params;
result.partition_type=partition_type;
result.resolution_range=resolution_range;

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

% def find_partition(N,sources,targets,
% weights=None, node_sizes=None, initial_membership=None, partition_type="RBC",
% resolution=1, n_iterations=-1, max_comm_size=0, rng_seed=None)

kwargs=pyargs('weights', weights, 'min_diff_bisect_value', params.min_diff_bisect_value, ...
    'min_diff_resolution', params.min_diff_resolution, 'number_iterations', int32(params.n_iterations), ...
    'linear_bisection',params.linear_bisection, 'rng_seed', int32(params.rng_seed));

res_profile=py.leiden.resolution_profile(N, sources, targets, ...
    partition_type, resolution_range, kwargs);

% result.res_profile=res_profile;
for i=1:length(res_profile)
    groups = py.numpy.array(res_profile{i}.membership);
    groups=double(groups);
    result.clusterID(:,i)=groups+1; %python is 0-based
    result.K(i)=max(groups+1); 
    result.mod(i)=res_profile{i}.modularity;
    result.q(i)=res_profile{i}.q;
    result.resolution(i)=res_profile{i}.resolution_parameter;
end


% groups = py.numpy.array(partition.membership);
% clusterID=double(groups);
% clusterID=clusterID+1; %python is 0-based
% K=max(clusterID); 

%relabel to be in decreasing size
% nT=arrayfun(@(x)nnz(clusterID==x),1:K);
% if params.doRelabel
%     [nT,ix]=sort(nT,'descend');
%     IDold=clusterID;
%     for i=1:K
%         clusterID(IDold==ix(i))=i;
%     end
% end
% result.K=K;
% result.clusterID=clusterID;
% result.clusterCounts = nT;