function merged_clust = merge_clusters(genes, counts, clust, options)
arguments
    genes
    counts
    clust
    options.simThr=0.5
    options.simMethod='overlap'
    options.nReps=1
    options.doPlot=0;
    options.cell_coords=[]
    options.tfidf_params=[]
end

for n=1:options.nReps
    
    cats=unique(clust);
    K=length(cats)
    
    if ~isempty(options.tfidf_params)
        TFIDF=gene_tfidf(genes,counts,clust,options.tfidf_params);
    else
        TFIDF=gene_tfidf(genes,counts,clust);
    end
    
    %fold-changes in expression, self vs notSelf
    dom=cell(1,K);
    for i=1:K
        dom{i}=TFIDF.("c"+string(i)).name;
    end
    
    %one-hot encode:
    uG=unique(cat(1,dom{:}));
    domG=false(K,length(uG));
    for i=1:K
        domG(i,ismember(uG,dom{i}))=true;
    end
    
    switch lower(options.simMethod)
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
    
    % merge and relabel
    [r,c]=find(Sim>options.simThr);
    merged_clust=clust;
    for i=1:length(r)
        tomerge=merged_clust==c(i);
        merged_clust(tomerge)=r(i);
    end
    
    remaincats=unique(merged_clust);
    newK=length(remaincats)
    nT=arrayfun(@(x)nnz(merged_clust==x),remaincats);
    [nT,ix]=sort(nT,'descend');
    IDold=merged_clust;
    for i=1:newK
        merged_clust(IDold==remaincats(ix(i)))=i;
    end
    
    if isequal(clust,merged_clust)
        break
    end
    
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
end