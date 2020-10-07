function [counts,gene_name,barcode,gene_id]=load10Xmtx(datapath, doSparse)
% load data from a folder containing 3 files (10X output)
% - count matrix (*matrix.mtx)
% - gene file (*genes.tsv)
% - barcode file (*barcodes.tsv)

barcodes_file=dir(datapath+"/*barcodes.tsv");
genes_file=dir(datapath+"/*genes.tsv");
if isempty(genes_file)
    genes_file=dir(datapath+"/*features.tsv");
end
matrix_file=dir(datapath+"/*matrix.mtx");

fid = fopen(fullfile(barcodes_file.folder,barcodes_file.name),'rt');
barcode = textscan(fid, '%s');
fclose(fid);
barcode=barcode{1};

fid = fopen(fullfile(genes_file.folder,genes_file.name),'rt');
geneData = textscan(fid, '%s %s');
fclose(fid);
gene_id=geneData{1};
gene_name=geneData{2};

% [mtxcounts,rows,colums,entries,representation,field,symmetry]...
%     =mmread(fullfile(matrix_file.folder,matrix_file.name));
counts=mmread(fullfile(matrix_file.folder,matrix_file.name));

if ~doSparse
    counts=full(counts);
end