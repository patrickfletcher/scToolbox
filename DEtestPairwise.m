function [p, adj_p, combs]=DEtestPairwise(X, group, method, keepTypes)
% just pairwise tests

correctionmethod='fdr';

minFrac=0;

nGenes=size(X,1);

group=categorical(group);
groupNames=categories(group);
groupNames=groupNames(:)';
groupCounts=countcats(group);
nGroups=length(groupNames);

if nGroups<2
    error('must have at least two groups to compare')
end

combs=nchoosek(groupNames,2);
nPairs=size(combs,1);

if exist('keepTypes')
    keepPair = any(ismember(combs,keepTypes),2);
    combs=combs(keepPair,:);
end

%filter genes if desired, else all genes tested
% keep=true(nGenes,1); %all
groupFrac=zeros(nGenes,nGroups);  %pass this in if available?
for i=1:nGroups
    groupFrac(:,i)=sum(X(:,group==groupNames{i})>0,2)/groupCounts(i);
end
keep=any(groupFrac > minFrac, 2); %any group expresses gene in at least some minimum fraction of cells (equiv to "either" in DEtest2)

keepix=find(keep);
nGenes2Test=nnz(keep);

% FDR is computed assuming all genes are present, to be conservative (?)
p=nan(nGenes,nPairs);

% tic
fprintf('\n')
tenth=round(nGenes2Test/10);
for i=1:nGenes2Test
    
    gix=keepix(i);
    for j=1:nPairs
        cells1=group==combs{j,1};
        cells2=group==combs{j,2};
        switch method
            case 'ttest'
                [~,this_p,this_ci,this_stats]=ttest2(X(gix,cells1),X(gix,cells2));
            case 'ranksum'
                [this_p,~,this_stats]=ranksum(X(gix,cells1),X(gix,cells2));
        end

        p(gix,j)=this_p;
    end
    if mod(i,tenth)==0
        fprintf('.')
    end
    
end
% toc

%give the smallest representable value to zero p-values so multiple comparisons has something to work with.
p(p==0)=eps(0);

%next the multiple comparisons per gene values
adj_p=ones(size(p));
for i=1:nPairs
    switch lower(correctionmethod)
        case 'bonferroni'
            %sequential bonferroni-holm
            adj_p(:,i)=bonf_holm(p(:,i));
        case {'fdr','bh'}
            %Benjamini-Hochberg FDR
            [~, ~, ~, adj_p(:,i)]=fdr_bh(p(:,i)); %method?
        otherwise
            adj_p=p;
    end
end

adj_p(adj_p>1)=1;