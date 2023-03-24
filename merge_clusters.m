function merged_clust = merge_clusters(genes, X, clust, options)
arguments
    genes
    X
    clust
    options.method='pca_distance' %degs

    %corr parameters
    options.summary_method="median" %mean
    options.distance_metric="correlation"

    %linkage for similarity clustering
    options.linkage="complete"

    %degs
    options.deg_method='median_cohen_d' %tfidf, median_delta_prop, ...
    options.blockvar=[]
    options.nTop = 20
    options.min_self_prop=0.1
    options.max_other_prop=1
    options.simMetric='overlap'

    %tfidf parameters
    options.tfidf_params=[]
    
    %general
    options.simThr=0.75
    options.nReps=1

    options.doPlot=0;
    options.plot_coords=[]
    options.plot_clusterSim=true
    options.plot_simThr=true;
end

thr = options.simThr;

%TODO:
% - support/force externally computed markers? what API makes sense?
% -- 1) do marker ID here. 2) pass markers in and just check similaritie here
% --- 1 offers iteration: can find markers of new clustering and do another merge
%
% - alternative cluster comparison methods: correlation of group means in
% PC space + linkage?
%
% - autothr?
clust0=clust;

for n=1:options.nReps

    cats = unique(clust);
    K=length(cats);
    cc = countcats(categorical(clust,cats));
    [~,ixs] = sort(cc);

    switch options.method
        case "degs"
            Sim = deg_sim(genes, X, clust, options);
        case "pca_distance"
            Sim = pca_distance(X, clust, metric=options.dist_metric);
    end

    y=squareform(Sim-eye(K),'tovector'); numel(y)
%     thr = geneThresholdOtsu(y(:)')

    % merge and relabel
    [r,c]=find(tril(Sim,-1)>=thr);

    merged_clust=clust;
%     uC=unique(c);
    newids=zeros(K,1);
    for i=1:K
        equiv = [i; r(c==i)];
        [~,mix]=max(cc(equiv));
        tomerge = ismember(merged_clust,cats(equiv));
        r(ismember(r,equiv))=equiv(mix);
        if length(equiv)>1 && any(tomerge)
            newids(equiv)=equiv(mix);
            merged_clust(tomerge)=cats(equiv(mix));
        end
        
%         allR=r(c==uC(i));
%         tomerge = ismember(merged_clust,cats(allR));
%         merged_clust(tomerge)=cats(uC(i));
    end

    remaincats=unique(merged_clust);
    newK=length(remaincats);
    oldK_newK=[K,newK]

    nT=arrayfun(@(x)nnz(merged_clust==x),remaincats);
    [nT,ix]=sort(nT,'descend');
    IDold=merged_clust;
    for i=1:newK
        merged_clust(IDold==remaincats(ix(i)))=i;
    end

    if isequal(clust,merged_clust)
        break
    end

    clust=merged_clust;
end


if options.doPlot
    fh=figure();clf
    ht=tiledlayout(1,3);

    oo=1:K;
    if options.plot_clusterSim
        Y=squareform(1-Sim);
        tree=linkage(Y,options.linkage);
        oo=optimalleaforder(tree,Y);
    end

    IM=Sim(oo,oo);

    cmap=turbo;
    if options.plot_simThr
        IM=IM-thr;
        cmap=split_cmap(Skip=16,nMid=1);
    end

    ax=nexttile(ht);
    imagesc(IM)
    colormap(ax,cmap);
    xticks(1:K)
    xticklabels(cats(oo))
    yticks(1:K)
    yticklabels(cats(oo))
    cb=colorbar;

    if options.plot_simThr
        CLIM=[-1,1]*max(abs(IM(:)));
        ax.CLim = CLIM;
        cb.Limits = [min(IM(:)),max(IM(:))];
    end

    ax=nexttile(ht);
    scatter_grp(options.plot_coords, clust0, textlabs=true, ax=ax);
    axis tight equal

    ax=nexttile(ht);
    scatter_grp(options.plot_coords, merged_clust, textlabs=true, ax=ax);
    axis tight equal
    drawnow
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sim, cats] = pca_distance(X, clustid, options)
arguments
    X 
    clustid
    options
end

[g,cats]=findgroups(clustid);
switch options.summary_method
    case "mean"
        clust_summary=splitapply(@mean,X,g(:));
    case "median"
        clust_summary=splitapply(@median,X,g(:));
end

Y=pdist(clust_summary, options.distance_metric);
D=squareform(Y);
Sim = 1-D./max(D(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sim, cats]= deg_sim(genes, expr, ident, options)
arguments
    genes
    expr
    ident
    options
end

switch options.deg_method
    case "tfidf"
        tfidf_params.N=options.nTop;
        tfidf_params.min_in_freq = options.min_self_prop;
        tfidf_params.max_out_freq = options.max_other_prop;
        tfidf_args=namedargs2cell(tfidf_params);
        [~,TFIDF]=gene_tfidf(genes,expr,ident,tfidf_args{:});
        M = table;
        M.celltype = TFIDF.clust_id;
        M.gene = TFIDF.name;

    case {'median_cohen_d','median_delta_prop'}
        summaryname=options.deg_method;
        effectname=strrep(summaryname,"median_","");
        ges=GroupedExpressionSummary(genes, expr, ident);
        ges.computeEffectSizes(effectname);
        ntop = options.nTop;
        minself = options.min_self_prop;
        maxother = options.max_other_prop;
        self=string(categories(ges.group));
        other=string(categories(ges.group));
        M=ges.topMarkers(summaryname,"top", self, other, nTop=ntop, ...
            min_self_prop=minself, max_other_prop=maxother);
end

cats = categories(categorical(M.celltype));
% [~,cats]=findgroups(M.celltype);
K = length(cats);

%one-hot encode gene lists:
uG=unique(M.gene);
isDEG=false(K,length(uG));
for i=1:K
    isDEG(i,ismember(uG,M.gene(M.celltype==cats{i})))=true;
end

switch lower(options.simMetric)
    case 'jaccard'
        simfun=@(A,B) nnz(A&B)/(nnz(A|B));
    case 'overlap'
        simfun=@(A,B) nnz(A&B)/min(nnz(A),nnz(B));  %asymmetric measure
    case 'dice'
        simfun=@(A,B) 2*nnz(A&B)/(nnz(A)+nnz(B));
    case 'correlation'
        simfun=@(A,B) corr(A(:),B(:));
    case 'nmi'
        simfun=@(A,B) nmi(A,B);
end

% now get the similarity matrix. Should this be symmetric? 

%pdist form: rows first, so can just do Sim(:) for linkage...
% Distances arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m – 1))

% asymmetric: overlap but only when row has smaller gene set?
Sim = zeros(K);
for j=1:K-1
    A=isDEG(j,:);
    for i=j+1:K
        B=isDEG(i,:);
        Sim(i,j)=simfun(A,B);
    end
end
% makes symmetric. could also use squareform(Sim(:))?
Sim = Sim + Sim' + eye(size(Sim));

end