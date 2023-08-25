function result = doPacMAP(X, knn, params, options)
arguments
    X
    knn = []
    params.nn_method = 'matlab'
    params.metric = "correlation"
    params.n_components = 2
    params.n_neighbors = 30
    params.num_iters = 450
    params.MN_ratio = 0.5
    params.FP_ratio = 2.0
    params.initY = "random"
    params.apply_pca = false
    params.rngSeed = 42
    options.verbose = false
    options.doPlot = false
end

n_obs = size(X,1);

% convert data for python
n_obs = int32(n_obs);
n_neighbors=int32(params.n_neighbors);
n_components=int32(params.n_components);
num_iters = int32(params.num_iters);
Xpy=py.numpy.array(X,'float32',pyargs('order','C'));
py_rand_state=int64(params.rngSeed);

%handle initY
initY=params.initY;
if isnumeric(initY) && numel(initY)<10
    %initY=[pc1,pc2] indicates index of PCs to use as initial points
    pcix=initY;
    signix=sign(pcix);
    pcix=abs(pcix);
    initY=X(:,pcix);
    initY(:,1)=signix(1)*initY(:,1);
    initY(:,2)=signix(2)*initY(:,2);
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));

elseif size(initY,1)==n_obs
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));

elseif initY~="spectral" && initY~="random" && initY~="pca"
    error('unknown initialization method for UMAP');
end

% core import
pacmap=py.importlib.import_module('pacmap');

% Matlab KNN (otherwise, let PacMAP do it)
if params.nn_method=="matlab"
    if ~isempty(knn) && size(knn.idx,1)==n_obs && params.n_neighbors<=size(knn.idx,2)
        knn_indices=knn.idx(:,1:params.n_neighbors);
        % knn_dists=knn.dists(:,1:params.n_neighbors);
    else
        disp('Computing neighbors (knnsearch)...')
        tic
        % k+1 because we don't want to remove self
        [knn_indices,~]=knnsearch(X,X,K=params.n_neighbors+1,Distance=metric); %quite fast!
        disp("knnsearch time: " + num2str(toc) + "s")
    end

    scaled_dist = py.numpy.ones(py.tuple([n_obs, n_neighbors])); % No scaling is needed;  
    scaled_dist = scaled_dist.astype(py.numpy.float32);
    nbrs = py.numpy.int32(knn_indices-1);
    pair_neighbors = pacmap.sample_neighbors_pair(Xpy, scaled_dist, nbrs, n_neighbors);

else
    pair_neighbors = string(missing); %None
end

% create PacMAP object and fit
disp('Computing PacMAP embedding...')
tic

embedding = pacmap.PaCMAP(n_components=n_components, n_neighbors=n_neighbors, num_iters=num_iters,...
    MN_ratio=params.MN_ratio, FP_ratio=params.FP_ratio, pair_neighbors=pair_neighbors, ...
    apply_pca=params.apply_pca, random_state=py_rand_state, distance=params.metric,...
    verbose=options.verbose);
embedded = embedding.fit_transform(Xpy, init=initY);

disp("PaCMAP time: " + num2str(toc) + "s")

result = params;
result.coords=double(embedded);   


% indptr = int64(emb_graph.indptr);
% r = int64(emb_graph.indices) + 1;
% data = double(emb_graph.data);
% 
% cdiff=diff(indptr);
% c=zeros(size(r), 'int64');
% for i=1:length(cdiff)
%     c(indptr(i)+1:indptr(i+1))=i*ones(cdiff(i),1, 'int64');
% end
% 
% result.graph=sparse(r, c, data, n_obs, n_obs);  %umap "connectivities"

% disp("extract and store result time: " + num2str(toc) + "s")

if options.doPlot
    fh=figure();clf
    scatter_grp(result.coords,ones(size(result.coords,1),1),gcols=0.66*[1,1,1], fig=fh);
    axis tight equal
    drawnow
end




%% %%%%%%%%%%%%%%%%%
%%%% PaCMAP:
% '''Pairwise Controlled Manifold Approximation.
% 
% Maps high-dimensional dataset to a low-dimensional embedding.
% This class inherits the sklearn BaseEstimator, and we tried our best to
% follow the sklearn api. For details of this method, please refer to our publication:
% https://www.jmlr.org/papers/volume22/20-1061/20-1061.pdf
% 
% Parameters
% ---------
% n_components: int, default=2
%     Dimensions of the embedded space. We recommend to use 2 or 3.
% 
% n_neighbors: int, default=10
%     Number of neighbors considered for nearest neighbor pairs for local structure preservation.
% 
% MN_ratio: float, default=0.5
%     Ratio of mid near pairs to nearest neighbor pairs (e.g. n_neighbors=10, MN_ratio=0.5 --> 5 Mid near pairs)
%     Mid near pairs are used for global structure preservation.
% 
% FP_ratio: float, default=2.0
%     Ratio of further pairs to nearest neighbor pairs (e.g. n_neighbors=10, FP_ratio=2 --> 20 Further pairs)
%     Further pairs are used for both local and global structure preservation.
% 
% pair_neighbors: numpy.ndarray, optional
%     Nearest neighbor pairs constructed from a previous run or from outside functions.
% 
% pair_MN: numpy.ndarray, optional
%     Mid near pairs constructed from a previous run or from outside functions.
% 
% pair_FP: numpy.ndarray, optional
%     Further pairs constructed from a previous run or from outside functions.
% 
% distance: string, default="euclidean"
%     Distance metric used for high-dimensional space. Allowed metrics include euclidean, manhattan, angular, hamming.
% 
% lr: float, default=1.0
%     Learning rate of the Adam optimizer for embedding.
% 
% num_iters: int, default=450
%     Number of iterations for the optimization of embedding. 
%     Due to the stage-based nature, we suggest this parameter to be greater than 250 for all three stages to be utilized.
% 
% verbose: bool, default=False
%     Whether to print additional information during initialization and fitting.
% 
% apply_pca: bool, default=True
%     Whether to apply PCA on the data before pair construction.
% 
% intermediate: bool, default=False
%     Whether to return intermediate state of the embedding during optimization.
%     If True, returns a series of embedding during different stages of optimization.
% 
% intermediate_snapshots: list[int], optional
%     The index of step where an intermediate snapshot of the embedding is taken.
%     If intermediate sets to True, the default value will be [0, 10, 30, 60, 100, 120, 140, 170, 200, 250, 300, 350, 450]
% 
% random_state: int, optional
%     Random state for the pacmap instance.
%     Setting random state is useful for repeatability.
% 
% save_tree: bool, default=False
%     Whether to save the annoy index tree after finding the nearest neighbor pairs.
%     Default to False for memory saving. Setting this option to True can make `transform()` method faster.
% '''

%%%% fit_transform:
% 
% '''Projects a high dimensional dataset into a low-dimensional embedding and return the embedding.
% 
% Parameters
% ---------
% X: numpy.ndarray
%     The high-dimensional dataset that is being projected. 
%     An embedding will get created based on parameters of the PaCMAP instance.
% 
% init: str, optional
%     One of ['pca', 'random']. Initialization of the embedding, default='pca'.
%     If 'pca', then the low dimensional embedding is initialized to the PCA mapped dataset.
%     If 'random', then the low dimensional embedding is initialized with a Gaussian distribution.
% 
% save_pairs: bool, optional
%     Whether to save the pairs that are sampled from the dataset. Useful for reproducing results.
% '''