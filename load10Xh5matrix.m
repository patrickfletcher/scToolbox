function [counts,genes,barcodes,geneIDs]=load10Xh5matrix(filename,h5path,doSparse, oldversion)
%load h5 compressed sparse column format into matlab full or sparse matrix

if ~exist('doSparse','var')
    doSparse=true;
end
if ~exist('oldversion','var')
    oldversion=false;
end

vals=double(h5read(filename,[h5path 'data']));
r=double(h5read(filename,[h5path 'indices']));
indptr=double(h5read(filename,[h5path 'indptr']));
shape=double(h5read(filename,[h5path 'shape']))';

r=r+1; %make sure it is 1 indexed for Matlab...

%column index must be expanded from indptr...
cdiff=diff(indptr);

% for building sparse matrix:
if doSparse
    c=zeros(size(r));
    for i=1:length(cdiff)
        c(indptr(i)+1:indptr(i+1))=i*ones(cdiff(i),1);
    end
    counts=sparse(r,c,vals,shape(1),shape(2));
else
    counts=zeros(shape);
    for i=1:length(cdiff)
        counts(r(indptr(i)+1:indptr(i+1)),i)=vals(indptr(i)+1:indptr(i+1));
    end
end

barcodes=h5read(filename,[h5path 'barcodes'])'; 
if oldversion
    geneIDs=h5read(filename,[h5path 'genes']);
    genes=h5read(filename,[h5path 'gene_names']);
else
    geneIDs=h5read(filename,[h5path 'features/id']);
    genes=h5read(filename,[h5path 'features/name']);
end

barcodes=string(barcodes(:));
genes=string(genes(:));
geneIDs=string(geneIDs(:));

%slight cleanup
genes=deblank(genes);
geneIDs=deblank(geneIDs);
barcodes=deblank(barcodes);

% for i=1:length(genes)
%     genes{i}=deblank(genes{i});
% end
% for i=1:length(geneIDs)
%     geneIDs{i}=deblank(geneIDs{i});
% end
% for i=1:length(barcodes)
%     barcodes{i}=deblank(barcodes{i});
% end