function [trimmed, info] = GOtrim(idtab, GO)
arguments
    idtab
    GO = []
end
% trim gene ontology list by traversing GO tree
% inspired by GO Trimming (https://doi.org/10.1186/1756-0500-4-267)
%
% input: table with columns "term_id", "intersections", "parents"
% returns a trimmed version of input

% minimal inputs: 
% - term id list, with GO term ids prefixed by "GO"
% - term intersections (list of genes per term)

% optionally:
% - ancestor information (is immediate parent enough? probably)
% * this must be only "is_a".

% gprof returns source, intersections, immediate parents.
% --> could already check immediate parents and trim via that...
% - PROBLEM: seems to return "regulates" closures in parents, which we
% don't want. We want only "is_a". Not even "part_of".

% find go terms
isGO=startsWith(idtab.term_id,'GO');
gotab=idtab(isGO,:);

uG=unique(cat(1,gotab.intersections{:}),'stable');

nTerms=height(gotab);
nUG=length(uG);

% need a binary term X gene matrix...
TGbin=false(nTerms,nUG);
for i=1:nTerms
    TGbin(i,:)=ismember(uG,gotab.intersections{i});
end

gsub=any(TGbin>0,1);
TGbin=TGbin(:,gsub);

% Similarity matrix between terms based on which genes are enriched
% - overlap coefficient using max(sizeA,sizeB) to normalize
% - create lower triangular matrix
X=zeros(nTerms,nTerms);
tic
for i=1:nTerms-1
    A=TGbin(i,:);
    nA=nnz(A);
    for j=i+1:nTerms  %this is slower
        B=TGbin(j,:);
        nB=nnz(B);
        nAB=nnz(A&B);
%         X(i,j)=nAB/min(nA,nB); %overlap of gene intersections
        X(j,i)=nAB/max(nA,nB); %overlap full matrix: assymetric
    end
end
% X=X+diag(ones(nTerms,1));

% if lower triangle value is 1, all larger term intersections were also
% present in smaller term. Edge case: term size is equal.
% L=tril(X,-1);
% [r,c]=find(L==1);

[r,c]=find(X==1);

largeT=gotab(r,:);
smallT=gotab(c,:);

if ~ismember('ancestors',smallT.Properties.VariableNames)
    smallT_ancestors=arrayfun(@(x)getancestors(GO,x,'Relationtype','is_a'),smallT.term_num,'UniformOutput',false);
else
    smallT_ancestors=smallT.ancestors;
end

%check if any "large" term is the parent of the corresponding "small" term
isParent=arrayfun(@(x,y)ismember(x,y{:}),largeT.term_num,smallT_ancestors);
% isParent=arrayfun(@(x,y)ismember(x,y{:}),largeT.term_id,smallT.parents); 

%don't consider self. only applies when results are from multi-queries
isSelf=largeT.term_num==smallT.term_num;
% isSelf=cellfun(@isequal, largeT.term_id, smallT.term_id);

% if all genes in larger are in smaller, larger term is a parent of smaller
% term, and the terms are not the same term, that means the larger term is
% redundant (less specific) and can be discarded.
removeix=isParent & ~isSelf;

remove=largeT(removeix,{'term_id','term_name','term_size'});
remove.child_id=smallT.term_id(removeix);
remove.child_name=smallT.term_name(removeix);
remove.child_size=smallT.term_size(removeix);

remove_ids=unique(remove.term_id);

redundant=ismember(idtab.term_id,remove_ids);
trimmed=idtab(~redundant,:);

info.isGO=isGO;
info.uG=uG;
info.TGbin=TGbin;
info.X=X;
info.largeT=largeT;
info.smallT=smallT;
info.isParent=isParent;
info.isSelf=isSelf;
info.removeix=removeix; %index into largeT/smallT
info.remove=remove;
info.remove_ids=remove_ids;
info.redundant=redundant;