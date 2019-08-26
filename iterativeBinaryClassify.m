function [classID,markers,markerScore,scoreThreshold,markerTable]=iterativeBinaryClassify(ncounts,tcounts,genes,initial_markers,params)
%iterative classifier given initial list of genes (of cells?)

%contraints/options: 
% - force keep initial gene list
% - max # genes
% - what if gene list shrinks to zero?

%initialize by gene list, or by cell subset?
% if cells: 1, DE => marker genes; 2, cells=classify; 3, repeat (1,2) until convergence
% if genes: 


GO=readtable('GO0000278_mitoticcellcycle.txt'); %too many...

DE_method='ranksum';
% DEmethod='kstest2';
DE_multcompare='fdr';
pthr=0.001;
genesToTest='group1';

% maxIter=10;
% selfThr=0.6;
% otherThr=0.05;
% keep_initial_ids=false;

% max_genes=100;
% %self/other comparison for dominant/specific
% % s.names={'self'}; %default
% s.minprct=20; %want to maximize this.
% s.fcPrctThr=2;
% s.fcExprThr=2;
% s.prctOrExpr=1;
% 
% % o.names={'other'};
% o.maxprct=10; %want to minimize this.  alt: sort by other% lowest to highest

%get the initial gene name indices
[gix,initial_markers]=getGeneIndices(initial_markers,genes.name); %remove any IDs not in our data set
markers=initial_markers; %could be list of string arrays, for multiclass support?
true_cells=false(size(tcounts,2));

%remove any genes not expressed in this cell set

iter=1;
while iter<=params.maxIter

    % classify based on new gene set

    geneThresholds=genes.thr(gix);
    genesAbove=tcounts(gix,:)>=geneThresholds;
    markerScore=sum(genesAbove,1);
    scoreThreshold=round(geneThresholdOtsu(markerScore));
%     scoreThreshold=length(gix)/3;
    new_true_cells=markerScore>=scoreThreshold;
    nnz(new_true_cells)
    
    if isequal(true_cells,new_true_cells) %check for stabilization of cyc cell set
        break
    else
        true_cells=new_true_cells;
    end
        
    %simple self/other method - advantage: doesn't constrain number of cells in self by imposing maxPrct in a cell type
    self=true_cells;
    other=~self;
    prct_self=sum(tcounts(:,self)>genes.thr,2)./nnz(self)*100;
    prct_other=sum(tcounts(:,other)>genes.thr,2)./nnz(other)*100;
%     mean_self=mean(ncounts(:,self),2);
%     mean_other=mean(ncounts(:,other),2);

    if params.doDE

        adj_p=DEtest2(DE_method,tcounts,self,other,DE_multcompare,genesToTest);
        keep=prct_self>params.selfPrctThr & prct_other<params.otherPrctThr & adj_p<pthr;

    else
        keep=prct_self>params.selfPrctThr & prct_other<params.otherPrctThr;
    end

%     keep=prct_self>params.selfPrctThr & prct_other<params.otherPrctThr & mean_self./mean_other>params.fcExprThr;

%     %dominantly expressed relative to identified cell types method
%     %first build a gene table with a column for "self", the identified cyc cells
%     self=true_cells;
%     prct_self=sum(tcounts(:,self)>genes.thr,2)./nnz(self)*100;
%     genetab=genes(:,contains(genes.Properties.VariableNames,'prct'));
%     genetab.prct_self=prct_self;
%     
% %     adj_p=DEtest2(DE_method,tcounts,self,other,DE_multcompare,genesToTest);
%     adj_p=ones(size(prct_self));
%     p.anova=adj_p;
%     p.mc=adj_p;
%     p.chi2=adj_p;
%     p.z=adj_p;
%     
%     [DOM,SPEC,dom,spec]=identifyDominantGenes(genetab,s,o,p);
    
    
    markers=genes.name(keep);
    markers=intersect(markers,GO.geneName);
    
    % constraints using known cell types? <thr% should be proliferating??
    
    
%     CYC=buildDEtable(self,other,genes,adj_p,N,T,genes,pthr,'up',1);
%     CYCspecific=filterDEtableFC(CYC,s,o);
%     marker_id=CYCspecific.gene;

%     if restrict_to_initial_ids
%         marker_id=intersect(marker_id,initial_id);
%     end
    if params.keep_initial_ids
        markers=unique([markers(:);initial_markers(:)]);
    end
    [gix,markers]=getGeneIndices(markers,genes.name);
    
    iter=iter+1;
end

markerTable=genes(gix,:);
classID=true_cells;