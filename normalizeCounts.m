function [A,sfs]=normalizeCounts(A,scale)
%normalize counts

%per cell
UMIPerGem=full(sum(A,1));

if ~exist('scale','var')||isempty(scale)
    scale=median(UMIPerGem);
end

% A=A./UMIPerGem.*scale; 
% A=A./repmat(UMIPerGem,length(A(:,1)),1).*scale;

%size factor (equivalent to above)
sfs=UMIPerGem/scale;

if issparse(A)
%     A=spfun(@(x)1./x,sfs);
    [i,j,a]=find(A);
    [m,n]=size(A);
    an=zeros(size(a));
    for k=1:length(a)
        an(k)=a(k)./sfs(j(k));
    end
    A=sparse(i,j,an,m,n);
else
    A=A./sfs;
end