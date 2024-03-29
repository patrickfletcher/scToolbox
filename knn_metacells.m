function MC = knn_metacells(C, scores, options)
arguments
    C %count matrix
    scores %reduced dimensionality scores (eg. PCA, MNN latent...)
    options.batch=[] %unused for now... do it externally
    options.precomputed=false
    options.n_neighbors=5
    options.metric='correlation'
    options.force_int_method='none'
    % options.method=@sum %must be simple functions for accumarray
end
% get KNN-aggregated counts for each cell ("meta-cells") - find knn of each
% cell, replace counts with the sum of counts across neighbors
% - make sure to include self.

% some sore of distance weighting? point is to use small K (sharing info in
% local neighborhood)

%TODO:  batch: find nearest neighbors only within within same batch. options:
% - subset scores to each batch. knn there, collate results.
% - for each batch, do pca + neighbors
% * no need if performed on integrated scores

if options.precomputed
    nnix=scores(:,1:options.n_neighbors); %assume self is removed??
else
    [nnix,~]=knnsearch(scores, scores, 'K',options.n_neighbors+1,'Distance',options.metric); %quite fast!
    nnix=nnix(:,2:end); %remove self
    % dist=dist(:,2:end);
end

%add counts of the K nearest neighbors to self
MC=C; %self
for i=1:options.n_neighbors
    MC=MC+C(:,nnix(:,i));
end


%force integer values? should be ints anyway
switch options.force_int_method
    case "round"
        MC = round(MC);
    case "ceil"
        MC = ceil(MC);
    otherwise
        %do nothing
end
