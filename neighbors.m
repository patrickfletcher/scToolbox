function result = neighbors(X, params, gparams, options)
arguments
    X
    params.n_neighbors=30
    params.metric="correlation"
    gparams.graph_type{mustBeMember(gparams.graph_type,["knn","snn"])}="knn"
    gparams.graph_weights = "adj"
    gparams.pruneThr = 0
    gparams.symmetrize = false
    options.makeGraph = false
    options.removeSelfNN = false
    options.verbose=false
end
%compute the k-nearest neighbors and their distances from a matrix of m
%from m observations in n dimensions, X [m rows by n columns]

%distance-weighted digraph object also stored: can get adjacency with
%builtin function

result = params;

disp('Computing neighbors...')
tic

[m,n]=size(X);
k=params.n_neighbors;

if ~options.removeSelfNN
    k = k+1;
end

% NOTE could save the searcher object for reuse if wanting to ajust
% n_neigbhors later

[idx,dists]=knnsearch(X,X,'K',k,'Distance',params.metric);

%Remove self distances
if options.removeSelfNN
    idx(:,1)=[];
    dists(:,1)=[];
end

disp("knnsearch time: " + num2str(toc) + "s")

%refinement options?
% - snn
% - snn + something to keep graph connected?
% https://github.com/LTLA/bluster/blob/master/src/build_snn.cpp

result.idx=idx;
result.dists=dists;

%graph representation
% TODO: split out graph building stuff into separate function?
if options.makeGraph
    sources = repelem(1:m, k)';
    targets = reshape(idx', [], 1);

    if gparams.graph_type=="snn" || gparams.graph_weights=="jaccard" || gparams.graph_weights=="number"
        [shared_neighbors, jaccard_index] = get_jaccard();
    end

    % - weight options? reciprocal knn_dist? "rank", "number", "jaccard" ???
    % https://rdrr.io/github/LTLA/bluster/man/makeSNNGraph.html
    switch gparams.graph_weights
        case "adj"
            weights = ones(size(sources));

        case "rank" 
            % w=k-r/2, where r is the smallest sum of ranks for any shared neighboring node
            % need to compute this alongside shared/jaccard... more complex cellfun function?

        case "number"
            weights = shared_neighbors;

        case "jaccard"
            weights = jaccard_index;

        case "dist" %reciprocal? 1/dist?  exponential kernel???
            weights = reshape(1./dists', 1, []);
            weights(isinf(weights))=0;
    end

    G = sparse(sources, targets, weights, m, m);

    switch gparams.graph_type
        case "knn"
            discards = false(m,1);

        case "snn" 
            discards = jaccard_index<gparams.pruneThr;
    end
    G(discards) = 0;

    %symmetrize? Louvain non-oriented does this anyway.
    if gparams.symmetrize
        G = (G + G')/2;
    end

    result.graph=G;
end

    %TODO: speed of looping? outer loop is on columns, cellfun loops on rows...
    function [shared_neighbors, jaccard_index] = get_jaccard()
        %similarity = Jaccard index (shared neighbors) -> what fraction of my neighbors include me as one of theirs?         
        shared_neighbors = zeros(size(idx));
        for i=1:k
            n_of_n=idx(idx(:,i),:);
            this_intersect = cellfun(@(c,d) length(intersect(c,d)), num2cell(idx, 2), num2cell(n_of_n, 2));
            shared_neighbors(:,i) = this_intersect;
        end
        shared_neighbors = reshape(shared_neighbors', [], 1);
        jaccard_index = shared_neighbors ./ (2*k - shared_neighbors); % <--- union = size(a) + size(b) - intersection
    end
end