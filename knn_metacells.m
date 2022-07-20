function MC = knn_metacells(C, scores, batch, params)
arguments
    C %count matrix
    scores %reduced dimensionality scores (eg. PCA, MNN latent...)
    batch = []
    params.precomputed=false
    params.n_neighbors=5
    params.metric='correlation'
end
% get KNN-aggregated counts for each cell ("meta-cells") - find knn of each
% cell, replace counts with the sum of counts across neighbors
% - make sure to include self.

% some sore of distance weighting? point is to use small K (sharing info in
% local neighborhood)

% batch: find nearest neighbors only within within same batch. options:
% - subset scores to each batch. knn there, collate results.
% - for each batch, do pca + neighbors
% * no need if performed on integrated scores

if params.precomputed
    nnix=scores(:,1:params.n_neighbors);
else
    nnix=knnsearch(scores,scores,'K',params.n_neighbors+1,'Distance',params.metric); %quite fast!
end

% accumarray?

MC=C;
for i=1:params.n_neighbors
    MC=MC+C(:,nnix(:,i));
end
