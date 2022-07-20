function result = neighbors(X, params, options)
arguments
    X
    params.n_neighbors=30
    params.metric="correlation"
    params.graph_type{mustBeMember(params.graph_type,["knn","snn"])}="knn"
    options.save_searcher
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

% NOTE could save the searcher object for reuse if wanting to ajust
% n_neigbhors later

% use k+1 because knnsearch returns self as 1st NN. 
[idx,dists]=knnsearch(X,X,'K',k+1,'Distance',params.metric);

%Remove self distances
idx(:,1)=[];
dists(:,1)=[];

disp("knnsearch time: " + num2str(toc) + "s")

%refinement options?
% - snn
% - snn + something to keep graph connected?

result.idx=idx;
result.dists=dists;

%for graph representation
% - weight options?

sources = repelem(1:m, k);
targets = reshape(idx', 1, []);
weights = reshape(dists', 1, []);

switch params.graph_type
    case "knn"
        G=digraph(sources,targets,weights); %already knn
    case "snn"
        %needs weights: "rank","number","jaccard"
        knnAdjacency = sparse(sources,targets,ones(size(sources)));
        A = knnAdjacency + knnAdjacency';
        A(A==1) = 0;
        % normalize to weights of 1 so the undirected mutual knn graph is unweighted
        A = A/2;
        G = graph(A);
end
result.graph=G;