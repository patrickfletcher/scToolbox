function [result, used_knn] = doUMAP(X, knn, params, options)
arguments
    X
    knn = [] 
    params.nn_method='matlab'
    params.metric="euclidean"
    params.n_neighbors=30
    params.n_components=2
    params.n_epochs=0  %let UMAP decide
    params.min_dist=0.3
    params.spread=1
    params.a=[]
    params.b=[]
    params.initial_alpha=1.0; %alpha is SGD learning rate
    params.gam=1.0; %gamma is repulsion strength
    params.negative_sample_rate=5; %double or int?
    params.initY = []
    params.init_noise = []
    params.init_offset = []
    params.init_norm = 'none'
    params.rngSeed=42  %default same seed... externally set if desired?
    params.do_densmap=false
    params.dens_lambda=2; %weight of density term in optimization
    params.dens_frac=0.3; %final fraction of epochs to use density term
    params.dens_var_shift=0.1; %small constant to prevent div by zero
    options.do_parallel=true;
    options.verbose=false
    options.doPlot = false
end

t0=tic;
n_obs = size(X,1); 

%TODO: is it correct to use the graph output by umap.fuzzy_simplicial_set
%for downstream clustering? --> the graph is what scanpy calls "connectivities"

%TODO - use of projection methods to add new points into UMAP?

% Cellranger defaults seem to be:
% npc=10, distance=correlation, n_neighbors=30, min_dist=0.3, seed=0.
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze
%
% Umap defaults to n_epochs=500 if n_obs<10000, else 200

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

isFullInitY = false;
if isempty(initY)
    initY="spectral";
    params.initY=initY;

elseif isnumeric(initY)&&numel(initY)<10
    isFullInitY = true;
    %initY=[pc1,pc2] indicates index of PCs to use as initial points
    pcix=initY;
    signix=sign(pcix);
    pcix=abs(pcix);
    initY=X(:,pcix);
    initY(:,1)=signix(1)*initY(:,1);
    initY(:,2)=signix(2)*initY(:,2);
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));

elseif size(initY,1)==n_obs
    isFullInitY = true;
    initY=py.numpy.array(initY,'float32',pyargs('order','C'));

elseif initY~="spectral" && initY~="random" && initY~="pca"
    error('unknown initialization method for UMAP');
end

if isFullInitY
    % scatter_grp(initY); axis on
    if ~isempty(params.init_offset)
        initY = initY + params.init_offset;
    end
    if ~isempty(params.init_noise)
        initY = initY + params.init_noise.*range(initY,2).*(rand(size(initY))-0.5);
    end
    switch params.init_norm
        case {'range','zscore','center','norm','scale','medianiqr'}
            initY = normalize(initY,1,params.init_norm);
        otherwise
            %nop
    end
    % scatter_grp(initY); axis on
end

if isstring(params.n_neighbors)||ischar(params.n_neighbors)
    switch params.n_neighbors
        case "sqrtN"
            params.n_neighbors=ceil(sqrt(n_obs));
        otherwise
            error('unknown n_neighbors string')
    end
end
if isempty(params.rngSeed)
    params.rngSeed = randi(intmax,1);
end
rng(params.rngSeed)


umap=py.importlib.import_module('umap.umap_');

%common python function arguments
n_neighbors=int64(params.n_neighbors);
metric_kwargs=py.dict();
py_rand_state=py.numpy.random.RandomState(int64(params.rngSeed));

% break into three steps for reuse of intermediates.
% 1. nearest neighbors --> nnIDX, nnDist (for impute!) - in python, quite slow! the bottleneck.

%matlab knn (otherwise, let UMAP do it)
if nn_method=="matlab"
    if ~isempty(knn) && size(knn.idx,1)==n_obs && params.n_neighbors<=size(knn.idx,2)
        doKNN=false;
        knn_indices=knn.idx(:,1:params.n_neighbors);
        knn_dists=knn.dists(:,1:params.n_neighbors);
    else
        disp('Computing neighbors (knnsearch)...')
        tic
        % k+1 because we don't want to remove self
        [knn_indices,knn_dists]=knnsearch(X,X,'K',params.n_neighbors+1,'Distance',metric); %quite fast!
        disp("knnsearch time: " + num2str(toc) + "s")
    end

    metric='precomputed'; %if I do KNN, UMAP doesn't use a metric
    kwargs = pyargs('knn_indices',py.numpy.int64(knn_indices-1),...
                    'knn_dists', py.numpy.array(knn_dists,'float32',pyargs('order','C')),...
                    'verbose',verbose,'return_dists',params.do_densmap);
else
    
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

    % disp('Computing neighbors (umap.nearest_neighbors)...')
    % kwargs = pyargs('angular',false,'verbose',true,'random_state',py_rand_state);
    % out=umap.nearest_neighbors(py.scipy.sparse.csr_matrix(X),n_neighbors,...
    %     metric,metric_kwargs,false,kwargs); 
    % knn_indices_py = out{1};
    % knn_dists_py=out{2};
    % knn_indices=double(knn_indices_py)+1; %+1 to convert to Matlab 1-based
    % knn_dists=double(knn_dists_py);
    % disp("nearest_neighbors time: " + num2str(toc) + "s")

    % metric='precomputed';
    % kwargs = pyargs('knn_indices',knn_indices_py,...
    %                 'knn_dists', knn_dists_py,...
    %                 'verbose',verbose,'return_dists',true);
    kwargs = pyargs('verbose',verbose,'return_dists',params.do_densmap);
end


disp('Computing fuzzy simplicial set...')

% 2. fuzzy simplicial set --> graph
tic
            
out_tuple=umap.fuzzy_simplicial_set(py.numpy.array(X,'float32',pyargs('order','C')),...
    n_neighbors, py_rand_state, metric, kwargs);
out_tuple=cell(out_tuple); %emb_graph, emb_sigmas, emb_rhos, emb_dists
emb_graph=out_tuple{1};
% emb_sigmas=out_tuple{2};
% emb_rhos=out_tuple{3}; %this looks like knn_dists to first neighbor (radius to first-NN)
% emb_dists=out_tuple{4}; %get this only if dens_map

disp("fuzzy_simplicial_set time: " + num2str(toc) + "s")

if ~isempty(params.a) && ~isempty(params.b)
    a=params.a;
    b=params.b;
else
    out=umap.find_ab_params(params.spread, params.min_dist);
    a=out{1}; 
    b=out{2};
end

if params.n_epochs >=0  
    % 3. embedding --> result
    disp('Computing simplicial set embedding...')
    tic
    
    n_components=int64(params.n_components);
    
    n_epochs=params.n_epochs;
    if n_epochs==0
        if n_obs<10000
            n_epochs=500;
        else
            n_epochs=200;
        end
    end
    n_epochs=int64(n_epochs);
    
    initial_alpha=params.initial_alpha; %alpha is SGD learning rate
    gam=params.gam; %gamma is repulsion strength
    negative_sample_rate=params.negative_sample_rate; %double or int?
    
    %do dens_map?
    densmap_kwds='';
    if params.do_densmap
        emb_dists=out_tuple{4};
        densmap_kwds=py.dict(pyargs('graph_dists',emb_dists,'lambda',params.dens_lambda,...
            'frac',params.dens_frac,'var_shift',params.dens_var_shift));
    end
    
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
    
    
    disp("simplicial set embedding time: " + num2str(toc) + "s")

end

%store the results
tic

result = params;
result.a=double(a);
result.b=double(b);
if params.n_epochs >= 0
    result.coords=double(embedded{1});   
end

% result.sigmas=double(emb_sigmas);
% result.rhos=double(emb_rhos); 
% result.dists=sparse(double(emb_dists.toarray()));

%just transfer the graph data, build sparse here. Much faster!
indptr = int64(emb_graph.indptr);
indices = int64(emb_graph.indices) + 1;
vals = double(emb_graph.data);
result.graph=build_sparse_matrix(vals, indptr, indices, n_obs, n_obs, 'csr');  %umap "connectivities"

if nargout==2
%connectivities used to generate UMAP's graph:
used_knn.indices=knn_indices;
used_knn.dists=knn_dists;
end

%connectivities in UMAP graph:
% [s,t,w]=find(cg.umap.graph);
% umap_graph_neighbors=splitapply(@(x) {x}, t, s);

disp("extract and store result time: " + num2str(toc) + "s")

if options.doPlot
    fh=figure();clf
    scatter_grp(result.coords,ones(size(result.coords,1),1),gcols=0.66*[1,1,1], fig=fh);
    axis tight equal
    drawnow
end


disp("doUMAP time: " + num2str(toc(t0)) + "s")