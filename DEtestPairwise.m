function [adj_pANOVA, adj_pMC, combs]=DEtestPairwise(X, group, method)
%this tests if all groups have same or different distributions (like 1-way anova with K-S two sample test)
% if not, multcompare is used to find out which are different

% method='kruskalwallis'; % could use anova, etc.

displayopt='off';
% ctype='tukey-kramer'; %has a floor of 10^-8
% ctype='dunn-sidak'; %has a floor, but also gives hard zeros
ctype='bonferroni'; %smooth values all the way to zero
% ctype='scheffe'; %smooth values all the way to zero
correctionmethod='fdr';

% if ~exist('minFrac','var')||isempty(minFrac)
    minFrac=0;
% end

nGenes=size(X,1);

group=categorical(group);
groupNames=categories(group);
groupCounts=countcats(group);
nGroups=length(groupNames);

if nGroups<2
    error('must have at least two groups to compare')
end

nPairs=nchoosek(nGroups,2);

%filter genes if desired, else all genes tested
% keep=true(nGenes,1); %all
groupFrac=zeros(nGenes,nGroups);  %pass this in if available?

for i=1:nGroups
    groupFrac(:,i)=sum(X(:,group==groupNames{i})>0,2)/groupCounts(i);
end
keep=any(groupFrac > minFrac, 2); %any group expresses gene in at least some minimum fraction of cells

keepix=find(keep);
nGenes2Test=nnz(keep);

% FDR is computed assuming all genes are present, to be conservative (?)
% should these be initialized to 1? I think so: expression=0 in all groups
pANOVA=ones(nGenes,1); 
pMC=ones(nGenes,nPairs); 
% C=nan(nGenes,6); %full comparisons matrix, includes estimate (difference between means?), CIs, p

tic
fprintf('\n')
tenth=round(nGenes2Test/10);
for i=1:nGenes2Test
    
    gix=keepix(i);
    
    switch method
        case 'anova1'
            [p,~,stats]=anova1(X(gix,:),group,displayopt);
        case 'kruskalwallis'
            [p,~,stats]=kruskalwallis(X(gix,:),group,displayopt);
    end
    pANOVA(gix,1)=p;

    %correction for nchoosek(nGroups,2) tests. 
    % 'alpha'??
    % 'Ctype': 'tukey-kramer' (default) | 'hsd' | 'lsd' | 'bonferroni' | 'dunn-sidak' | 'scheffe'
    c=multcompare(stats,'display',displayopt,'Ctype',ctype);

    pMC(gix,:)=c(:,end)';
    
    if mod(i,tenth)==0
        fprintf('.')
    end
    
end
toc

groupNames=groupNames(:)'; %force row vector to get simple pair compare right:
combs=groupNames(c(:,1:2));

%give the smallest representable value to zero p-values so multiple comparisons has something to work with.
pANOVA(pANOVA==0)=eps(0);
pMC(pMC==0)=eps(0);

%correction for nGenes tests
%first, the ANOVA p values
switch lower(correctionmethod)
    case 'bonferroni'
        %sequential bonferroni-holm
        adj_pANOVA=bonf_holm(pANOVA);
    case {'fdr','bh'}
        %Benjamini-Hochberg FDR
        [~, ~, ~, adj_pANOVA]=fdr_bh(pANOVA); %method?
    otherwise
        adj_pANOVA=pANOVA;
end

adj_pANOVA(adj_pANOVA>1)=1;

%next the multiple comparisons per gene values
adj_pMC=ones(size(pMC));
for i=1:nPairs
    switch lower(correctionmethod)
        case 'bonferroni'
            %sequential bonferroni-holm
            adj_pMC(:,i)=bonf_holm(pMC(:,i));
        case {'fdr','bh'}
            %Benjamini-Hochberg FDR
            [~, ~, ~, adj_pMC(:,i)]=fdr_bh(pMC(:,i)); %method?
        otherwise
            adj_pMC=pMC;
    end
end

adj_pMC(adj_pMC>1)=1;