function [counts,genes,barcodes, fileInfo]=load10Xh5matrix(filename, doSparse, isPreV3)
arguments
    filename
    doSparse=0
    isPreV3=0
end
%load h5 compressed sparse column format into matlab full or sparse matrix

hi = h5info(filename); h5path = hi.Groups(:).Name; h5path(end+1)='/';

if ~exist('doSparse','var'), doSparse=true; end
if ~exist('preVersion3','var'), isPreV3=false; end

vals=int32(h5read(filename,[h5path 'data']));
r=int64(h5read(filename,[h5path 'indices']));
indptr=int64(h5read(filename,[h5path 'indptr']));
shape=double(h5read(filename,[h5path 'shape']))';

r=r+1; %make sure it is 1 indexed for Matlab...

%column index must be expanded from indptr...
cdiff=diff(indptr);

% for building sparse matrix:
if doSparse
    c=zeros(size(r),'like',r);
    for i=1:length(cdiff)
        c(indptr(i)+1:indptr(i+1))=i*ones(cdiff(i),1);
    end
    counts=sparse(r,c,vals,shape(1),shape(2));
else
    counts=zeros(shape,'like',vals);
    for i=1:length(cdiff)
        counts(r(indptr(i)+1:indptr(i+1)),i)=vals(indptr(i)+1:indptr(i+1));
    end
end

if nargout>=2
    genes = table;
    if isPreV3
        gene_id=h5read(filename,[h5path 'genes']);
        gene_name=h5read(filename,[h5path 'gene_names']);
    else
        gene_id=h5read(filename,[h5path 'features/id']);
        gene_name=h5read(filename,[h5path 'features/name']);
        gene_feature_type=h5read(filename,[h5path 'features/feature_type']);
        gene_genome=h5read(filename,[h5path 'features/genome']);
    end
    gene_id=string(gene_id(:));
    genes.id=deblank(gene_id);
    gene_name=string(gene_name(:));
    genes.name=deblank(gene_name);
    
    if ~isPreV3
        gene_feature_type=string(gene_feature_type(:));
        genes.feature_type=deblank(gene_feature_type);
        gene_genome=string(gene_genome(:));
        genes.genome=deblank(gene_genome);
    end
end

if nargout>=3
    barcodes=h5read(filename,[h5path 'barcodes'])'; 
    barcodes=string(barcodes(:));
    barcodes=deblank(barcodes);
end

if nargout==4
    fileInfo = hi; 
end