function [adj_pANOVA, adj_pMC, combs]=DEtestPairwiseAnova(X, group, options)
arguments
    X
    group
    options.Method='kruskalwallis'
    options.anova_correction='bonferroni'
    options.P_correction='fdr'
    options.MinFrac=0
    options.MinCells=0
end
%this tests if all groups have same or different distributions (like 1-way anova with K-S two sample test)
% if not, multcompare is used to find out which are different

% options.anova_correction
%  'tukey-kramer'; %has a floor of 10^-8
%  'dunn-sidak'; %has a floor, but also gives hard zeros
%  'bonferroni'; %smooth values all the way to zero
%  'scheffe'; %smooth values all the way to zero

%NOTE: if 2 groups, pANOVA is same as 2-sample test (eg. kruskalwallis
%pANOVA = p from ranksum). Should fall back to simply doing the latter.

nGenes=size(X,1);

group=categorical(group);
groupNames=categories(group);
groupCounts=countcats(group);
nGroups=length(groupNames);

if nGroups<2
    error('must have at least two groups to compare')
end

groupNames=groupNames(:)';
combs=nchoosek(groupNames,2);
nPairs=size(combs,1);

%filter genes if desired
groupCells=zeros(nGenes,nGroups);
for i=1:nGroups
    groupCells(:,i)=sum(X(:,group==groupNames{i})>0,2);
end
groupFrac=groupCells./groupCounts';
keep=any(groupCells > options.MinCells & groupFrac > options.MinFrac, 2); %any group expresses gene in at least some minimum fraction of cells (equiv to "either" in DEtest2)

keepix=find(keep);
nGenes2Test=nnz(keep);

% FDR is computed assuming all genes are present, to be conservative (?)
% should these be initialized to 1 or nan?
pANOVA=nan(nGenes,1); 
pMC=nan(nGenes,nPairs); 
% mean_ranks=zeros(nGenes,nPairs);
% C=nan(nGenes,6); %full comparisons matrix, includes estimate (difference between means?), CIs, p

% tic
fprintf('testing %d genes', nGenes2Test)
tenth=floor(nGenes2Test/10);
for i=1:nGenes2Test
    gix=keepix(i);
    switch options.Method
        case 'anova1'
            [p,~,stats]=anova1(X(gix,:),group,'off');
%         case 'wanova'
%             [p,F,df1,df2]=wanova(X(gix,:),group);
%             stats.
            %Games-Howell post hoc test ?
        case 'kruskalwallis'
            [p,~,stats]=kruskalwallis(X(gix,:),group,'off');
        otherwise
            error('unsupported test method')
    end
    pANOVA(gix,1)=p;

    %correction for nchoosek(nGroups,2) tests. 
    % 'Ctype': 'tukey-kramer' (default) | 'hsd' | 'lsd' | 'bonferroni' | 'dunn-sidak' | 'scheffe'
    if nGroups>2
        c=multcompare(stats,'Ctype',options.anova_correction,'Display','off');
        pMC(gix,:)=c(:,end)';
    end
    
    if mod(i,tenth)==0
        fprintf('.')
    end
    
end
% toc
if nGroups==2
    pMC=pANOVA;
end

%give the smallest representable value to zero p-values so multiple comparisons has something to work with.
pANOVA(pANOVA==0)=eps(0);
pMC(pMC==0)=eps(0);

%correction for nGenes tests
%first, the ANOVA p values
switch lower(options.P_correction)
    case 'bonferroni'
        %sequential bonferroni-holm
        adj_pANOVA=bonf_holm(pANOVA);
    case {'fdr','bh'}
        %Benjamini-Hochberg FDR
        [~, ~, ~, adj_pANOVA]=fdr_bh(pANOVA); %method?
    otherwise
        adj_pANOVA=pANOVA;
end


%next the multiple comparisons per gene values
adj_pMC=ones(size(pMC));
if nGroups>2
    for i=1:nPairs
        switch lower(options.P_correction)
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
end

adj_pANOVA(adj_pANOVA>1)=1;
adj_pMC(adj_pMC>1)=1;

%for simple two-sample test, kruskal-wallace p = ranksum p. The
%multiple-comparison by multcompare doesn't make sense?
if nGroups==2
    adj_pMC=adj_pANOVA;
end