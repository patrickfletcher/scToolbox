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

disp('Performing Leidenalg clustering...')

%add path to python script if needed...
Pypath = py.sys.path;
MLpath=string(path).split(';');
sctoolpath=MLpath(contains(MLpath,'scToolbox'));
if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
    insert(Pypath,int32(0),sctoolpath);
end

pyleiden = py.importlib.import_module('leiden');

%X should contain the graph from umap... need to convert to a csr matrix
%can't pass sparse to python

if contains(class(X),"graph")
    X=adjacency(X);
end

N=size(X,1);
[sources,targets,weights]=find(X);
N = int32(N); %cast to int for python
sources = int32(sources-1); %convert to zero-based index as well
targets = int32(targets-1);

n_part =length(params.resolution);

counts={};
K=zeros(1,n_part);
mod=zeros(1,n_part);
q=zeros(1,n_part);
clusterID=zeros(size(X,1),n_part);
for r=1:n_part

    % def find_partition(N,sources,targets,
    % weights=None, node_sizes=None, initial_membership=None, partition_type="RBC",
    % resolution=1, n_iterations=-1, max_comm_size=0, rng_seed=None)
    
    kwargs=pyargs('weights', weights, 'resolution', params.resolution(r),'n_iterations', ...
        int32(params.n_iterations), 'max_comm_size', int32(params.max_comm_size), ...
        'rng_seed', int32(params.rng_seed));
    
    pypart=pyleiden.find_partition(N,sources,targets,partition_type, kwargs);
    
    groups = py.numpy.array(pypart.membership);
    ids=double(groups);
    ids=ids+1; %python is 0-based
    k=max(ids); 
    
    %relabel to be in decreasing size
    nT=arrayfun(@(x)nnz(ids==x),1:k);
    if params.doRelabel
        [nT,ix]=sort(nT,'descend');
        IDold=ids;
        for i=1:k
            ids(IDold==ix(i))=i;
        end
    end

    K(r)=k;
    mod(r)=pypart.modularity;
    q(r)=pypart.q;
    clusterID(:,r)=ids;
    counts{r} = nT;
end

result.K=K;
result.mod=mod;
result.q=q;
result.clusterID=clusterID;
result.counts=counts;

pinfo=table;
pinfo.res=params.resolution(:);
pinfo.K=K(:);
pinfo.mod=mod(:);
pinfo.q=q(:);

result.info=pinfo;