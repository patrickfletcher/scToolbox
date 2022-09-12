function [tree, D, oo] = build_cluster_tree(X, clustid, options)
arguments
    X
    clustid
    options.summary_method="median"
    options.distance_metric="correlation"
    options.linkage="complete"
    options.reorder_factors=true
    options.do_plot=false
end
%use data in X - expression or PC scores - to create a linkage tree between
%clusters specified in clustid. Computes distance matrix D from "average"
%representative of each cluster, in X-space.

% [g,ID]=findgroups(clustid);
[g,ID]=grp2idx(string(clustid));

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

S = 1-D./max(D(:));

if options.do_plot
%     figure()
%     imagesc(S(oo,oo))
%     set(gca,"Ydir","normal")
%     xticks(1:length(ID));
%     yticks(1:length(ID));
%     xticklabels(ID(oo))
%     yticklabels(ID(oo))
%     axis(gca,"square")
    hcg=clustergram(S, Symmetric=0, Colormap=turbo, ColumnLabels=ID, RowLabels=ID,...
        linkage=options.linkage, ColumnPDist=options.distance_metric, RowPDist=options.distance_metric);
%     axis(gca,"square")
end
        