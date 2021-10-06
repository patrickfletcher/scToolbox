function score = score_genes(query_genes, expr, gene_pool, options)
arguments
    query_genes
    expr
    gene_pool
    options.ctrl_size = 50
    options.n_bins = 25
    options.min_cells = 1  %TODO - this is unused
    options.rng_seed = 42
end
% compute gene scores as done by AddModuleScore in Seurat and
% scanpy.tl.score_genes: Tirosh2016 DOI: 10.1126/science.aad0501

rng(options.rng_seed)

% cellsExpr=sum(expr>0,2);
% expr=expr(cellsExpr>=min_cells,:);
% gene_pool=gene_pool(cellsExpr>=min_cells);

meanExpr=mean(expr,2);
[sortedMean,ixs]=sort(meanExpr); %genebin_ids should index into sorted ixs

%compute background scores in each expression bin
bin_edges=round(linspace(1,length(gene_pool),options.n_bins+1));
[genebin_id, bin_edges]=discretize(1:numel(gene_pool),bin_edges);

%compute query gene scores
qix=getGeneIndices(query_genes, gene_pool);
qixs=getGeneIndices(query_genes, gene_pool(ixs));
qbins=unique(genebin_id(qixs));

bg_ix=[];
for i=1:length(qbins)
    this_bin=find(genebin_id==qbins(i));
    rand_ix=randsample(this_bin,options.ctrl_size); 
%     rand_ix=this_bin(randperm(length(this_bin),options.ctrl_size));          %could use randsample instead.
%     rand_ix=this_bin(randi([1,length(this_bin)],1,options.ctrl_size));
    bg_ix=[bg_ix;ixs(rand_ix)]; %gene IDs from the this bin
end
bg_ix=setdiff(bg_ix,qix); %don't include query genes

q_means=mean(expr(qix,:),1); %mean of expression across query genes, per cell
bg_means=mean(expr(bg_ix,:),1);

score=q_means - bg_means;