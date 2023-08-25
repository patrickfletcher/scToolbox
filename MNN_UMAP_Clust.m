function cg = MNN_UMAP_Clust(datafile, csub, splitby, genes, cells, mnnopts, opts)
arguments
    datafile
    csub
    splitby {string,char,cellstr}
    genes
    cells
    mnnopts.gene_subset=false
    mnnopts.min_mean=0.1 
    mnnopts.do_pooledsizefactors=false
    mnnopts.do_multibatch=true
    mnnopts.min_mean_hvg=0.1 
    mnnopts.do_densityweights=false
    mnnopts.n_features=0
    mnnopts.fdr_thr=1
    mnnopts.var_thr=0
    mnnopts.mergeo=[]
    mnnopts.d=50
    mnnopts.k=20
    mnnopts.prop_k=[]
    mnnopts.ndist=3
    opts.metric="correlation"
    opts.knn=30
    opts.n_neighbors=30
    opts.n_epochs=0
    opts.min_dist=0.3
    opts.initY="spectral"
    opts.do_3D=true
    opts.initY3d="spectral"
    opts.resolution=1
    opts.do_plot=false
end

if islogical(csub)
    cg.subset=csub(:);
else
    cg = csub; %it's an existing scd - overwrite the fields
end 

csubtab=cells(:,["id", mnnopts.splitby]); 
csubtab.keep=cg.subset; %ix into cleaned data
% datafile, cellsub, splitby, genes
cg.mnn=scran_fastMNN(datafile, csubtab, splitby, k=mnnopts.k, d=mnnopts.d, ...
    do_pooledsizefactors=mnnopts.do_pooledsizefactors, do_multibatch=mnnopts.do_multibatch, ...
    do_densityweights=mnnopts.do_densityweights, fdr_thr=mnnopts.fdr_thr, ...
    n_features=mnnopts.n_features,merge_order=mnnopts.mergeo, verbose=1);

cg.mknn=neighbors(cg.mnn.coords, n_neighbors=opts.knn, metric=opts.metric, removeSelfNN=false);

cg.umap=doUMAP(cg.mnn.coords, cg.mknn, n_neighbors=opts.n_neighbors, ...
    n_epochs=opts.n_epochs, metric=opts.metric, min_dist=opts.min_dist, initY=opts.initY);
if opts.do_3D
    cg.umap3=doUMAP(cg.mnn.coords, cg.mknn, n_components=3, n_neighbors=opts.n_neighbors, ...
        n_epochs=opts.n_epochs, metric=opts.metric, min_dist=opts.min_dist, initY=opts.initY3d);
end

cg.clust=leiden(cg.umap.graph, "RBC", resolution=opts.resolution);

if opts.do_plot
    scatter_grp(cg.umap.coords, cg.clust.clusterID(:,1), textlabs=true);
    axis tight equal
end