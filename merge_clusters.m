function merged_clust = merge_clusters(genes, counts, clust, options)
arguments
    genes
    counts
    clust
%     options.deg_method='tfidf'
    options.deg_params=[]
    options.simThr=0.5
    options.simMetric='overlap'
%     options.nReps=1
    options.doPlot=0;
    options.cell_coords=[]
end

%TODO:
% - support/force externally computed markers? what API makes sense?
% -- 1) do marker ID here. 2) pass markers in and just check similaritie here 
% --- 1 offers iteration: can find markers of new clustering and do another merge 
%
% - alternative cluster comparison methods: correlation of group means in
% PC space + linkage?

% for n=1:options.nReps
    
    cats=unique(clust);
    K=length(cats);
    
    %self vs notSelf DEGs
    if ~isempty(options.deg_params)
        deg_par_args=namedargs2cell(options.deg_params);
        DEG=gene_tfidf(genes,counts,clust,deg_par_args{:});
    else
        DEG=gene_tfidf(genes,counts,clust);
    end
    
    %extract lists of genes
    deg=cell(1,K);
    for i=1:K
        deg{i}=DEG.("c"+string(i)).name;
    end
%     cellfun(@length,deg)
    
    %one-hot encode gene lists:
    uG=unique(cat(1,deg{:}));
    domG=false(K,length(uG));
    for i=1:K
        domG(i,ismember(uG,deg{i}))=true;
    end
    
    switch lower(options.simMetric)
        case 'jaccard'
            simfun=@(A,B) nnz(A&B)/(nnz(A|B));
        case 'overlap'
            simfun=@(A,B) nnz(A&B)/min(nnz(A),nnz(B));
        case 'dice'
            simfun=@(A,B) 2*nnz(A&B)/(nnz(A)+nnz(B));
        case 'correlation'
            simfun=@(A,B) corr(A(:),B(:));
        case 'nmi'
            simfun=@(A,B) nmi(A,B);
    end
            
    % now get the similarity matrix
    Sim = zeros(K);
    for i=1:K-1
        A=domG(i,:);
        for j=i+1:K
            B=domG(j,:);
            Sim(i,j)=simfun(A,B);
        end
    end

%     Sim=geneset_similarities(M.gene,M.celltype,options.simMetric);

    %options for other merge method? linkage? 
    
    % merge and relabel
    [r,c]=find(Sim>options.simThr);
    merged_clust=clust;
    for i=1:length(r)
        tomerge=merged_clust==c(i);
        merged_clust(tomerge)=r(i);
    end
    
    %TODO>use rename_clusters.
    remaincats=unique(merged_clust);
    newK=length(remaincats);
    oldK_newK=[K,newK]
    nT=arrayfun(@(x)nnz(merged_clust==x),remaincats);
    [nT,ix]=sort(nT,'descend');
    IDold=merged_clust;
    for i=1:newK
        merged_clust(IDold==remaincats(ix(i)))=i;
    end
    
%     if isequal(clust,merged_clust)
%         break
%     end
    
    if options.doPlot
        fh=figure(1);clf
        ht=tiledlayout(1,3);
        
        ax=nexttile(ht);
        imagesc(Sim)
        xticks(1:K)
        yticks(1:K)
        colorbar
        
        ax=nexttile(ht);
        plotScatter(options.cell_coords,'group',clust,[],fh,ax,[],'index');
        axis tight equal
        
        ax=nexttile(ht);
        plotScatter(options.cell_coords,'group',merged_clust,[],fh,ax,[],'index');
        axis tight equal
        drawnow
    end
    
    clust=merged_clust;
% end