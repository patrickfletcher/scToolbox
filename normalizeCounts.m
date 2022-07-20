function [ncounts,sfs]=normalizeCounts(counts, options)
arguments
    counts
    options.method="libsize"
    options.scale=[]
    options.max_frac=[]
end
%normalize counts

%TODO: center size factors around 1??

%TODO: other methods???
% - what if cell borrowed info from k-nn? needs way to get KNN (pca???)
% -- 1. naive counts_per_cell, naive scale --> sfs0=counts/median(counts)
%  [knn: sfs0 -> norm, logcounts -> pca -> knn.idx]
% -- 2. summed counts_per_cell: countsK, summed scale=median(countsK): scaleK --> sfsK=countsK/scaleK 
% -- 3. sfs --> sfsK*

%TODO: params.method, etc. 

method=options.method;
scale=options.scale;
max_frac=options.max_frac;

exclude_hiexp=~isempty(max_frac);

%sum counts per cell
counts_per_cell=full(sum(counts,1));

%scale computed before exclude genes
if isempty(scale)
    scale=median(counts_per_cell);
end

if exclude_hiexp
    %redo the counts_per_cell
    high_frac=full(sum(counts>max_frac*counts_per_cell,2));
    counts_per_cell=full(sum(counts(~high_frac,:),1));
    scale=median(counts_per_cell);
end

%size factors
sfs=counts_per_cell/scale;

% %alt: DESeq - do I need to geomean only positive vals to avoid divbyzeros?
% gm=geomean(counts,2); %pseudoreference sample (cell)
% sfs=median(counts./gm,1);

if issparse(counts)
%     A=spfun(@(x)1./x,sfs);
    [i,j,a]=find(counts);
    [m,n]=size(counts);
    an=zeros(size(a));
    for k=1:length(a)
        an(k)=a(k)./sfs(j(k));
    end
    ncounts=sparse(i,j,an,m,n);
else
    ncounts=counts./sfs;
end