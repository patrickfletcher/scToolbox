function [counts,genes,barcodes,geneIDs]=load10XData(name,doSparse)
% load data from a folder containing 3 files (10X output)
% - count matrix (name_matrix.mtx)
% - gene file (name_genes.tsv)
% - barcode file (name_barcodes.tsv)
%

%TODO: merge the h5 and MTX versions together into one function
fid = fopen([name,'_genes.tsv'],'rt');
geneData = textscan(fid, '%s %s');
fclose(fid);

geneIDs=geneData{1};
genes=geneData{2};

fid = fopen([name,'_barcodes.tsv'],'rt');
barcodes = textscan(fid, '%s');
fclose(fid);
barcodes=barcodes{1};

[counts,rows,colums,entries,representation,field,symmetry]=mmread([name,'_matrix.mtx']);

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