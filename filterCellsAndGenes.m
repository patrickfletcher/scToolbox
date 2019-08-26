function [counts,genes,barcodes,geneIDs,genesubset]=filterCellsAndGenes(counts,genes,barcodes,geneIDs,cellsubset,minCells)

%minimum inputs: counts, genes, cellsubset

if ~exist('minCells','var')||~isempty(minCells)
    minCells=1;
end

% subset cells
counts=counts(:,cellsubset);

% subset genes
cellsWithGene=sum(counts>0,2); %count how many cells expressed each gene
genesubset=cellsWithGene>=minCells;

counts=counts(genesubset,:);
genes=genes(genesubset);

if exist('barcodes','var')&&~isempty(barcodes)
    barcodes=barcodes(cellsubset);
else
    barcodes=[];
end

if exist('geneIDs','var')&&~isempty(geneIDs)
    geneIDs=geneIDs(genesubset);
else
    geneIDs=[];
end


