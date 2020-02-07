function M=knnGraph(X,k,type,pruneThr,symmetrize)
% similarity matrix using knn

%TODO: pruning?

% [idx,Dist]=knnsearch(X,X,'k',k);
[idx,Dist]=knnsearch(X,X,'k',k+1);
idx(:,1) = []; Dist(:,1) = []; %remove self-distances

if ~exist('type','var')
    type='adjacency';
end

if ~exist('pruneThr','var')
    pruneThr=1/15;
end

if ~exist('symmetrize','var')
    symmetrize=false;
end

% Create the adjacency matrix for linkage
M = zeros(size(X,1));
for ii = 1:length(X)
    
    %simple Knn
    switch lower(type)
        case 'adjacency'
            %adjacency matrix only
            M(ii,idx(ii,:)) = 1;
            
        case 'euclidean'
            %similarity = 1/Euclidean distance
%             M(ii,idx(ii,:)) = Dist(ii,:); 
            M(ii,idx(ii,:)) = 1./Dist(ii,:); 
            
        case 'jaccard'
            %similarity = Jaccard index (shared neighbors)
            pt_neighbs = idx(ii,:);
            n_of_n = idx(pt_neighbs,:);
            shared_neighbors = sum(ismember(n_of_n, pt_neighbs),2);
            % intersection and union sum to k
            weights = shared_neighbors ./ (2*k-shared_neighbors); % Jaccard coefficient
            %pruning? seems like if weight<pruneKNN, set weight=0.
            weights(weights<pruneThr)=0;
            M(ii,idx(ii,:)) = weights;
    end
end

%symmetrize? Louvain non-oriented does this anyway.
if symmetrize
    M=M+M';
end

% D=full(knn2jaccard(idx)); %lower triangular, D+D' is equal to my symmetrized.