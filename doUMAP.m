function result = doUMAP(X, params, knn, figID, valuenames, colors)

% Cellranger defaults seem to be:
% npc=10, distance=correlation, n_neighbors=30, min_dist=0.3, seed=0.
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze
%
% Differences: I use Matlab knnsearch, euclidean
%
% Umap defaults to n_epochs=500 if size(X,1)<10000, else 200

disp('Computing UMAP...')

rng(params.rngSeed)

initY=params.initY;
if isempty(initY)
    initY="spectral";
    params.initY=initY;
elseif isnumeric(initY)&&numel(initY)==2
    %initY=[pc1,pc2] indicates index of PCs to use as initial points
    pcix=initY;
    signix=sign(pcix);
    pcix=abs(pcix);
    initY=X(:,pcix);
    initY(:,1)=signix(1)*initY(:,1);
    initY(:,2)=signix(2)*initY(:,2);
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));
elseif size(initY,1)==size(X,1)
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));
elseif initY~="spectral" && initY~="random"
    error('unknown initialization method for UMAP');
end

if ~isnumeric(params.n_neighbors) && params.n_neighbors=="sqrtN"
    params.n_neighbors=floor(sqrt(size(X,1)));
end

%initialize results
result.params = params;

doKNN=true;
if exist('knn','var') && ~isempty(knn) && params.n_neighbors<=size(knn.indices,2)
    doKNN=false;
    knn_indices=knn.indices(:,1:params.n_neighbors);
    knn_dists=knn.dists(:,1:params.n_neighbors);
end
    
%common python function arguments
n_neighbors=int64(params.n_neighbors);
verbose=true;
metric='euclidean';
% metric='correlation';
% metric='cosine';
metric_kwargs=py.dict();
py_rand_state=py.numpy.random.RandomState(int64(params.rngSeed));

% break into three steps for reuse of intermediates.

% 1. nearest neighbors --> nnIDX, nnDist (for impute!) - in python, quite slow! the bottleneck.
%matlab knn
if doKNN
    tic
    [knn_indices,knn_dists]=knnsearch(X,X,'K',params.n_neighbors,'Distance',metric); %quite fast!
    % knn_indices(:,1)=[]; %remove self distances
    % knn_dists(:,1)=[];
    disp("knnsearch time: " + num2str(toc) + "s")
%     result.knn_indices=knn_indices;
%     result.knn_dists=knn_dists;
end

%python version
% do_sparse=1; %seems much faster for py nearest neighbors
% if do_sparse %transform seemed to fail when X is sparse.
%     X=py.scipy.sparse.csr_matrix(X);
% end
% tic
% kwargs = pyargs('verbose',true,'random_state',py_rand_state);
% out=py.umap.umap_.nearest_neighbors(py.scipy.sparse.csr_matrix(X),n_neighbors,...
%     metric,metric_kwargs,false,kwargs); 
% result.knn_indices=double(out{1})+1; %+1 to convert to Matlab 1-based
% result.knn_dists=double(out{2});
% toc


% 2. fuzzy simplicial set --> graph
tic
kwargs = pyargs('knn_indices',py.numpy.int64(knn_indices-1),...
                'knn_dists', py.numpy.array(knn_dists,'float32',pyargs('order','C')),...
                'verbose',verbose);
out_tuple=py.umap.umap_.fuzzy_simplicial_set(py.numpy.array(X,'float32',pyargs('order','C')),...
    n_neighbors, py_rand_state, metric, kwargs);
out_tuple=cell(out_tuple);
fuzzy_simplicial_set=out_tuple{1};
result.graph=sparse(double(fuzzy_simplicial_set.toarray())); %is this correct?
disp("fuzzy_simplicial_set time: " + num2str(toc) + "s")


% 3. embedding --> result
tic
n_components=int64(2);
initial_alpha=1.0; %alpha is SGD learning rate
gam=1.0; %gamma is repulsion strength
negative_sample_rate=5; %double or int?
n_epochs=int64(params.n_epochs);
init=initY;

out=py.umap.umap_.find_ab_params(params.spread, params.min_dist);
a=out{1}; b=out{2};
result.a=double(a);
result.b=double(b);

kwargs = pyargs('verbose',verbose);
embedded = py.umap.umap_.simplicial_set_embedding(X, fuzzy_simplicial_set,...
    n_components, initial_alpha, a, b, gam, negative_sample_rate,... 
    n_epochs, init, py_rand_state, metric, metric_kwargs, kwargs);
result.coords=double(embedded);
disp("simplicial_set_embedding time: " + num2str(toc) + "s")
    
% %just run umap directly
% kwargs = pyargs('n_neighbors',int32(params.n_neighbors),...
%     'min_dist',params.min_dist,'spread',params.spread,...
%     'n_epochs',int32(params.n_epochs),...
%     'n_components',int32(2),'verbose',verbose,'init',initY,...
%     'random_state',int32(params.rngSeed));
% 
% umap_obj=py.umap.UMAP(kwargs).fit(X);
% % res=umap_obj.transform(X); %this is for transforming new data.
% result.coords=double(umap_obj.embedding_);
% result.graph=sparse(double(umap_obj.graph_.toarray())); %just save the scipy.csr_matrix to be passed to leiden

if exist('figID','var')
    figure(figID);clf
    plotScatter(result.coords,'value',valuenames,colors,figID);
    axis tight
    axis equal
    drawnow
end
