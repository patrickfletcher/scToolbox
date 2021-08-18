function [p, adj_p, combs]=DEtestPairwise(X, group, options)
arguments
    X
    group
    options.Method='ranksum'
    options.MinFrac=0
    options.MinCells=0
    options.Groups=[]
%     options.Genes=[]
end
% just pairwise tests

correctionmethod='fdr';

[nGenes,nCells]=size(X);

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

if ~isempty(options.Groups)
    keepPairs = any(ismember(combs,options.groups),2);
    combs=combs(keepPairs,:);
end

cells1=false(nPairs,nCells);
cells2=false(nPairs,nCells);
for j=1:length(nPairs)
    cells1(j,:)=group==combs{j,1};
    cells2(j,:)=group==combs{j,2};
end

%filter genes
groupCells=zeros(nGenes,nGroups);
for i=1:nGroups
    groupCells(:,i)=sum(X(:,group==groupNames{i})>0,2);
end
groupFrac=groupCells./groupCounts';
keep=any(groupCells > options.MinCells & groupFrac > options.MinFrac, 2); %any group expresses gene in at least some minimum fraction of cells (equiv to "either" in DEtest2)

keepix=find(keep);
nGenes2Test=nnz(keep);

% FDR is computed assuming all genes are present, to be conservative (?)
p=nan(nGenes,nPairs);

% tic
fprintf('testing %d genes', nGenes2Test)
tenth=floor(nGenes2Test/10);
for i=1:nGenes2Test
    gix=keepix(i);
    for j=1:nPairs
        switch options.Method
            case 'ttest'
                [~,this_p]=ttest2(X(gix,cells1(j,:)), X(gix,cells2(j,:)));
            case 'ranksum'
                [this_p]=ranksum(X(gix,cells1(j,:)), X(gix,cells2(j,:)));
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