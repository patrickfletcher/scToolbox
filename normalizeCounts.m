function [ncounts,sfs]=normalizeCounts(counts, scale, max_frac)
%normalize counts

%TODO: params.method, etc. 

exclude_hiexp=false;
if exist('max_frac','var')
    exclude_hiexp=true;
end

%sum counts per cell
counts_per_cell=full(sum(counts,1));

%scale computed before exclude genes
if ~exist('scale','var')||isempty(scale)
    scale=median(counts_per_cell,1);
end

if exclude_hiexp
    %redo the counts_per_cell
    high_frac=full(sum(counts>max_frac*counts_per_cell,2));
    counts_per_cell=full(sum(counts(~high_frac,:),1));
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