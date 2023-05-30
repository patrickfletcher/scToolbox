function result = mixing_score(coords, group, K, nnix, nndists, options)
arguments
    coords
    group
    K = 100
    nnix = []
    nndists = []
    options.method = "entropy"
    options.weigth_by_dist=false
    options.knnmetric="correlation"
end

% find K-nearest neighbors if needed
n_knn = size(nnix,2);
if K>n_knn
    [nnix,nndists]=knnsearch(coords,coords,'K',K+1,'Distance',options.knnmetric);
    nnix(:,1)=[];
    nndists(:,1)=[];
elseif K<n_knn    
    nnix(:,K+1:end)=[];
    nndists(:,K+1:end)=[];
end
group=removecats(group(:)); 
groups=categories(group);
nGroups=length(groups);

%initialize the results structure with the parameters
result=options;
result.K=K;

%proportions expected by chance
P0 = countcats(group)'./size(group,1);
result.P0=P0;

%all methods based on proportion of K-NN of each batch
nn_groups=group(nnix);
g_cnts=countcats(nn_groups,2);
Pi=g_cnts./K;
result.Pi=Pi;

Wi=ones(size(Pi));
if options.weigth_by_dist
    % weighting should favor closer neighbors. (following cms score)
    % sum(weights) should = 1...
    % -- could weight each observation, or use mean (as cms does)
    % compute mean distance to nearest neighbors of same type.

    % adjust P0 by overall mean nndist:
    W0=mean(nndists,2);
    P0=P0.*W0;

    for i=1:nGroups
        nnd=nndists; 
        nnd(nn_groups~=groups{i})=nan; %keep this groups' distances only
        Wi(:,i)=mean(1./nnd,2,'omitnan'); %can end up NaN if no NN are this group
    end
    Wi(isnan(Wi))=0; %adjust the nan weights. the Pi for those is zero anyway...
    Wi=Wi./sum(Wi,2);
    result.Wi=Wi;
    Pi=Pi.*Wi;
end

% other options? 
% - proportion tests? Chi2: (O-E)2^2/E
% - distribution tests (e.g. cms nn-distance distributions)? hyge?
switch options.method
    case "max_dev"
        %1 - sum(absolute deviations from P0)
%         dev=Pi-P0;
        dev=abs(Pi-P0)./P0;
        scores = max(dev,[],2);
        scores = 1 - scores/max(scores);
%         result.dev=dev;

    case "sum_dev"
        %1 - sum(absolute deviations from P0)
%         dev=Pi-P0;
        dev=abs(Pi-P0)./P0;
        scores = sum(dev,2);
        scores = 1 - scores/max(scores);
%         result.dev=dev;

    case "chi2"
        chi2 = sum((Pi-P0).^2./P0,2);
        dof=nGroups-1;
%         scores = 1-chi2/max(chi2);
        p = chi2cdf(chi2,dof,'upper');
        scores = p;
%         scores = chi2;
%         result.dev=dev;

    case "entropy"
        % normalized Shannon entropy of the set of proportions
        maxE=log(nGroups);
        PlogP=Pi.*log(Pi); 
        PlogP(Pi==0)=0; %replace undefined value with function's limit
        scores = -1/maxE*sum(PlogP,2);

    case "LISI"
        LISI = 1./sum(Pi.^2,2);
%         scores=LISI;
        scores = (LISI-1)/(nGroups-1);

end

result.scores=scores;