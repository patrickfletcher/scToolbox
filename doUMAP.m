function [result, used_knn] = doUMAP(X, knn, params, options)
arguments
    X
    knn = [] 
    params.nn_method='matlab'
    params.metric="correlation"
    params.n_neighbors=30
    params.n_components=2
    params.n_epochs=0  %let UMAP decide
    params.min_dist=0.3
    params.spread=1
    params.initial_alpha=1.0; %alpha is SGD learning rate
    params.gam=1.0; %gamma is repulsion strength
    params.negative_sample_rate=5; %double or int?
    params.initY = []
    params.rngSeed=42
    params.do_densmap=false
    params.dens_lambda=2; %weight of density term in optimization
    params.dens_frac=0.3; %final fraction of epochs to use density term
    params.dens_var_shift=0.1; %small constant to prevent div by zero
    options.do_parallel=false;
    options.verbose=false
    options.figID = []
end

t0=tic;

%TODO: is it correct to use the graph output by umap.fuzzy_simplicial_set
%for downstream clustering?

%TODO - use of projection methods to add new points into UMAP?

% Cellranger defaults seem to be:
% npc=10, distance=correlation, n_neighbors=30, min_dist=0.3, seed=0.
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze
%
% Umap defaults to n_epochs=500 if size(X,1)<10000, else 200

%TODO: might be better to use the highest level interface, UMAP(), to
%prevent frequent breaking of code..
%
% kwargs = pyargs('n_neighbors',int32(params.n_neighbors),...
%     'min_dist',params.min_dist,'spread',params.spread,...
%     'n_epochs',int32(params.n_epochs),...
%     'n_components',int32(2),'verbose',verbose,'init',initY,...
%     'random_state',int32(params.rngSeed));
% 
% umap_obj=umap.UMAP(kwargs).fit(X);
% % res=umap_obj.transform(X); %this is for transforming new data.
% result.coords=double(umap_obj.embedding_);
% result.graph=sparse(double(umap_obj.graph_.toarray())); %just save the scipy.csr_matrix to be passed to leiden

%for now unpack a few params:
initY=params.initY;
metric=params.metric;
nn_method=params.nn_method;
verbose=options.verbose;
do_parallel=options.do_parallel;

if isempty(initY)
    initY="spectral";
    params.initY=initY;
elseif isnumeric(initY)&&numel(initY)<10
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
elseif initY~="spectral" && initY~="random" && initY~="pca"
    error('unknown initialization method for UMAP');
end

if isstring(params.n_neighbors)||ischar(params.n_neighbors)
    switch params.n_neighbors
        case "sqrtN"
            params.n_neighbors=round(sqrt(size(X,1)));
        otherwise
            error('unknoqn n_neighbors string')
    end
end

%initialize results (flat struct with all params + results)
result = params;

rng(params.rngSeed)

%common python function arguments
n_neighbors=int64(params.n_neighbors);

% break into three steps for reuse of intermediates.
% 1. nearest neighbors --> nnIDX, nnDist (for impute!) - in python, quite slow! the bottleneck.
doKNN=true;
if ~isempty(knn) && size(knn.idx,1)==size(X,1) && params.n_neighbors<=size(knn.idx,2)
    doKNN=false;
    knn_indices=knn.idx(:,1:params.n_neighbors);
    knn_dists=knn.dists(:,1:params.n_neighbors);
end

%matlab knn
if doKNN && nn_method=="matlab"
    disp('Computing neighbors...')
    tic
    [knn_indices,knn_dists]=knnsearch(X,X,'K',params.n_neighbors+1,'Distance',metric); %quite fast!
%     knn_indices(:,1)=[]; %remove self distances
%     knn_dists(:,1)=[];
    disp("knnsearch time: " + num2str(toc) + "s")
    metric='precomputed'; %if I do KNN, UMAP doesn't use a metric
end

metric_kwargs=py.dict();
py_rand_state=py.numpy.random.RandomState(int64(params.rngSeed));

%python version
% do_sparse=1; %seems much faster for py nearest neighbors
% if do_sparse %transform seemed to fail when X is sparse.
%     X=py.scipy.sparse.csr_matrix(X);
% end
% tic
% kwargs = pyargs('verbose',true,'random_state',py_rand_state);
% out=umap.nearest_neighbors(py.scipy.sparse.csr_matrix(X),n_neighbors,...
%     metric,metric_kwargs,false,kwargs); 
% result.knn_indices=double(out{1})+1; %+1 to convert to Matlab 1-based
% result.knn_dists=double(out{2});
% toc


umap=py.importlib.import_module('umap.umap_');
disp('Computing UMAP...')

% 2. fuzzy simplicial set --> graph
tic
if nn_method=="matlab" || ~doKNN
kwargs = pyargs('knn_indices',py.numpy.int64(knn_indices-1),...
                'knn_dists', py.numpy.array(knn_dists,'float32',pyargs('order','C')),...
                'verbose',verbose,'return_dists',true);
else
kwargs = pyargs('verbose',verbose,'return_dists',true);
end
            
out_tuple=umap.fuzzy_simplicial_set(py.numpy.array(X,'float32',pyargs('order','C')),...
    n_neighbors, py_rand_state, metric, kwargs);
out_tuple=cell(out_tuple); %emb_graph, emb_sigmas, emb_rhos, emb_dists
emb_graph=out_tuple{1};
% emb_sigmas=out_tuple{2};
% emb_rhos=out_tuple{3}; %this looks like knn_dists to first neighbor (radius to first-NN)
% emb_dists=out_tuple{4}; %get this only if dens_map

disp("fuzzy_simplicial_set time: " + num2str(toc) + "s")


% 3. embedding --> result
tic

n_components=int64(params.n_components);

n_epochs=params.n_epochs;
if n_epochs==0
    if size(X,1)<10000
        n_epochs=500;
    else
        n_epochs=200;
    end
end
n_epochs=int64(n_epochs);

initial_alpha=params.initial_alpha; %alpha is SGD learning rate
gam=params.gam; %gamma is repulsion strength
negative_sample_rate=params.negative_sample_rate; %double or int?

out=umap.find_ab_params(params.spread, params.min_dist);
a=out{1}; b=out{2};

%do dens_map?
densmap_kwds='';
if params.do_densmap
    emb_dists=out_tuple{4};
    densmap_kwds=py.dict(pyargs('graph_dists',emb_dists,'lambda',params.dens_lambda,...
        'frac',params.dens_frac,'var_shift',params.dens_var_shift));
end

%doesn't work on my system for some reason
% do_parallel=true;
if do_parallel
    py_rand_state=py.numpy.random.RandomState();
end


kwargs = pyargs('verbose',verbose,'densmap',params.do_densmap,'densmap_kwds', ...
    densmap_kwds,'output_dens',false, 'parallel',do_parallel);

%metric here is for Spectral initY only
metric=params.metric; %set it back to the original metric choice in case we used pre-computed.
% - has special code branch for Euclidean - is it faster??
embedded = umap.simplicial_set_embedding(X, emb_graph,...
    n_components, initial_alpha, a, b, gam, negative_sample_rate,... 
    n_epochs, initY, py_rand_state, metric, metric_kwargs, kwargs);


disp("simplicial_set_embedding time: " + num2str(toc) + "s")
    
%store the results
tic
result.graph=sparse(double(emb_graph.toarray())); %is this correct?
% result.sigmas=double(emb_sigmas);
% result.rhos=double(emb_rhos); 
% result.dists=sparse(double(emb_dists.toarray())); 
result.a=double(a);
result.b=double(b);
result.coords=double(embedded{1});
disp("extract and store result time: " + num2str(toc) + "s")

if nargout==2
%connectivities used to generate UMAP's graph:
used_knn.indices=knn_indices;
used_knn.dists=knn_dists;
end

%connectivities in UMAP graph:
% [s,t,w]=find(cg.umap.graph);
% umap_graph_neighbors=splitapply(@(x) {x}, t, s);

disp("doUMAP time: " + num2str(toc(t0)) + "s")

if ~isempty(options.figID)
    figure(options.figID);clf
    scatter_grp(result.coords,ones(size(result.coords,1),1),gcols=0.66*[1,1,1],fig=options.figID);
    axis tight equal
    drawnow
end
