function [tree, D, oo] = build_cluster_tree(X, clustid, options)
arguments
    X
    clustid
    options.summary_method="mean"
    options.distance_metric="euclidean"
    options.linkage="ward"
    options.reorder_factors=true
end
%use data in X - expression or PC scores - to create a linkage tree between
%clusters specified in clustid. Computes distance matrix D from "average"
%representative of each cluster, in X-space.

[g,ID]=findgroups(clustid);

switch options.summary_method
    case "mean"
        clust_summary=splitapply(@mean,X,g(:));
    case "median"
        clust_summary=splitapply(@mean,X,g(:));
end

Y=pdist(clust_summary,options.distance_metric);
D=squareform(Y);
tree=linkage(Y,options.linkage);
oo=optimalleaforder(tree,Y);

% imagesc(D(oo,oo))
        