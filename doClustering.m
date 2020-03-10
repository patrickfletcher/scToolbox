function result = doClustering(scores, params)
% tic
result = params;
switch params.method
    case 'kmeans'
        disp('Performing K-means clustering...')
        clusterID=kmeans(scores,params.K);
        
    case 'linkage'
        disp('Performing hierarchical clustering...')
        Z = linkage(scores,params.linkage,params.metric);
        clusterID = cluster(Z,'Maxclust',params.K);
        result.Z=Z; %the tree
        
    case 'modularity'
        disp('Performing graph-based clustering...')
        %Seurat uses KNN=10 and pruning=1/15 as defaults...
        D=knnGraph(scores,params.KNN,params.metric, params.pruneThr);
        COMTY = cluster_jl_cpp(D,1,1,0,0);
        clusterID=COMTY.COM{end};
        result.K=length(COMTY.SIZE{end});
end

%relabel to be in decreasing size
nT=arrayfun(@(x)nnz(clusterID==x),1:result.K);
[nT,ix]=sort(nT,'descend');
IDold=clusterID;
for i=1:result.K
    clusterID(IDold==ix(i))=i;
end
result.clusterID=clusterID;
result.clusterCounts = nT;
% clusttime=toc