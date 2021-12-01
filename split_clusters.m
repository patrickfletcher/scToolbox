function new_clust=split_clusters(coords, clust, tosplit, options)
arguments
    coords
    clust
    tosplit
    options.splitmethod='umap'
    options.clustmethod
    options.resolution
    options.n_neighbors
    options.min_dist
    options.n_epochs
end

ix2split=ismember(clust, tosplit);
sub_umap=doUMAP(coords(ix2split,:), umapopts);

sub_clust=cellgroup.clust.(rval);
sub_clust.resolution=0.5;
sub_clust=doClustering(sub_umap.graph,sub_clust);

splitclust(ix2split)=sub_clust.clusterID+cellgroup.clust.(rval).K;

new_K=length(unique(splitclust));

%relabel to be in decreasing size
nT=arrayfun(@(x)nnz(splitclust==x),1:max(splitclust));
[nT,ix]=sort(nT,'descend');
IDold=splitclust;
for i=1:new_K
    splitclust(IDold==ix(i))=i;
end

new_clust=cellgroup.clust.(rval);
new_clust.K=new_K;
new_clust.clusterID=splitclust;
new_clust.clusterCounts=nT;
new_clust.split=tosplit;
new_clust.sub_umap=sub_umap;
new_clust.sub_clust=sub_clust;