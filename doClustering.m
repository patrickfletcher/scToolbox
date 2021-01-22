function result = doClustering(scores, params)
% tic
result = params;
switch params.method
    case 'kmeans'
        disp('Performing K-means clustering...')
        clusterID=kmeans(scores,params.K);
        result.K=params.K;
        
    case 'linkage'
        disp('Performing hierarchical clustering...')
        Z = linkage(scores,params.linkage,params.metric);
        clusterID = cluster(Z,'Maxclust',params.K);
        result.Z=Z; %the tree
        result.K=params.K; 
        
    case 'modularity'
        disp('Performing graph-based clustering...')
        %Seurat uses kNN=10 and pruning=1/15 as defaults...
        D=knnGraph(scores,params.kNN,params.metric, params.pruneThr);
        COMTY = cluster_jl_cpp(D,1,1,0,0);
        clusterID=COMTY.COM{end};
        result.K=length(COMTY.SIZE{end});
    case 'leiden'
        disp('Performing Leidenalg (python) clustering...')
        
        %add path to python script if needed...
        Pypath = py.sys.path;
        MLpath=string(path).split(';');
        sctoolpath=MLpath(contains(MLpath,'scToolbox'));
        if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
            insert(Pypath,int32(0),sctoolpath);
        end
        
        %scores should contain the graph from umap... need to convert to a csr matrix
%         graph = full(scores); %can't pass sparse to python
        N=size(scores,1);
        [sources,targets,weights]=find(scores);
        N = int32(N); %cast to int for python
        sources = int32(sources-1); %convert to zero-based index as well
        targets = int32(targets-1);
        clusterID=double(py.leiden.run(N,sources,targets,weights,params.resolution));
%         clusterID=double(py.leiden.run(graph,params.resolution));
        clusterID=clusterID+1; %python is 0-based
        result.K=max(clusterID); 
end

if ~isfield(params, 'doRelabel')
    params.doRelabel=true;
end

%relabel to be in decreasing size
nT=arrayfun(@(x)nnz(clusterID==x),1:result.K);
if params.doRelabel
    [nT,ix]=sort(nT,'descend');
    IDold=clusterID;
    for i=1:result.K
        clusterID(IDold==ix(i))=i;
    end
end
result.clusterID=clusterID;
result.clusterCounts = nT;
% clusttime=toc