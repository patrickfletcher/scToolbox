function p_adj=DEtest2(method, X, group1, group2, correctionmethod, genesToTest, minFrac)
%two sample differential expression test with multiple corrections
% tests each row, comparing two groups of columns
%method specifies which test to use

%group1 = logical, group2 absent/empty => group2=~group1
%group1 = logical, group2=logical => two sample test

if ~exist('method','var')||isempty(method)
    method='ranksum';
end

if ~exist('group2','var')||isempty(group2)
    group2=~group1;
end

if ~exist('correctionmethod','var')||isempty(correctionmethod)
    correctionmethod='fdr';
end

if ~exist('minFrac','var')||isempty(minFrac)
    minFrac=0;
end

if ~exist('genesToTest','var')||isempty(genesToTest)
    genesToTest='either';
end

[nGenes,nCells]=size(X);

if length(group1)~=nCells
    error('grouping variable must contain one element per column of X');
end

%filter genes if desired, else all genes tested 
frac1=sum(X(:,group1)>0,2)/nnz(group1); 
frac2=sum(X(:,group2)>0,2)/nnz(group1);
switch genesToTest
    case 'all'
        keep=true(nGenes,1); %all
    case 'group1'
        keep=frac1 > minFrac;
    case 'group2'
        keep=frac2 > minFrac;
    case 'either'
        keep=frac1 > minFrac | frac2 > minFrac;
end

keepix=find(keep);
nGenes2Test=nnz(keep);

% FDR is computed assuming all genes are present, to be conservative (?)
p=nan(nGenes,1); 

tic
% fprintf('\n')
tenth=ceil(nGenes2Test/10);
for i=1:nGenes2Test

    gix=keepix(i);
    switch method
        case 'ranksum'
            p(gix,1)=ranksum(X(gix,group1),X(gix,group2));
        case 'kstest2'
            [~,p(gix,1)]=kstest2(X(gix,group1),X(gix,group2));
        case 'ttest2'
            [~,p(gix,1)]=ttest2(X(gix,group1),X(gix,group2));
    end

    if mod(i,tenth)==0 || i==nGenes2Test
        fprintf('.')
    end
    
end

%give the smallest representable value to zero p-values so multiple comparisons has something to work with.
p(p==0)=eps(0);
toc

%correction for nGenes tests
switch lower(correctionmethod)
    case 'bonferroni'
        %sequential bonferroni-holm
        p_adj=bonf_holm(p);
    case {'fdr','bh'}
        %Benjamini-Hochberg FDR
        [~, ~, ~, p_adj]=fdr_bh(p); %method?
end

p_adj(p_adj>1)=1;