function result=doDimRedAndCluster(scdata,ncounts,tcounts,genes,params_hvg,params_pca,params_tsne, params_clust, doPlot, group, colors)
% HVGs, PCA, tSNE, clustering
% result is a struct with fields hvg, pca, tsne, and clust

%TODO: Cellranger uses all genes but truncates to fixed # PCs. Here, I use
%fixed number HVG, measure # PCs by perm test... Which is correct?

if ~isempty(scdata) && isstruct(scdata)
    result=scdata; %appends results to existing struct
    ncounts=ncounts(scdata.genesub,scdata.subset);
    tcounts=tcounts(scdata.genesub,scdata.subset);
    gene_name=genes.name(scdata.genesub);
    gene_id=genes.id(scdata.genesub);
else
    gene_name=genes.name;
    gene_id=genes.id;
end

%store the parameters used in the same struct as results
hvg_results=params_hvg;
pca_results=params_pca;
tsne_results=params_tsne;
clust_results=params_clust;

if ~exist('group','var')||isempty(group)
    group=ones(1,size(ncounts,2));
end
if ~exist('colors','var')||isempty(colors)
    colors=[0.5,0.5,0.5];
end
if ~exist('doPlot','var')||isempty(doPlot)
    doPlot=false;
end

% 1. hvgs
tic
disp('Finding highly variable genes...')
% params_hvg.nHVG=min(round(size(ncounts,2)/2),params_hvg.nHVG);

if doPlot
    figure(1);clf
    hvgix=findVariableGenes(ncounts,params_hvg,1);
    drawnow
else
    hvgix=findVariableGenes(ncounts,params_hvg);
end
hvg_results.ix=hvgix;
hvg_results.gene_name=gene_name(hvgix);
hvg_results.gene_id=gene_id(hvgix);

hvgtime=toc;

% 2. pca
tic

X=tcounts(hvgix,:);

%is rescaling necessary? covariance vs correlation? logtransformed counts are roughly on same scale, in [0,3]
switch params_pca.scale_method
    case 'zscore'
        % Seurat scaling
        X=(X-mean(X,2))./std(X,[],2); %zscore gene distributions
        X(X>params_pca.maxScaled)=params_pca.maxScaled;
        X(X<-params_pca.maxScaled)=-params_pca.maxScaled;  %
        
    case 'unit'
        %unitize, then subtract mean
        X=(X-min(X,[],2))./(max(X,[],2)-min(X,[],2));
        X=(X-mean(X,2));
        
    case 'center'
        X=(X-mean(X,2));
        
    case 'none'
end

if params_pca.npc<1
    % Find number of significant PCs
    disp('Performing PCA permutation test...')
    p=findSignificantPCs(X',params_pca.permute_reps,0.05);
    nPCs=find(~(p==1/params_pca.permute_reps),1,'first')-1;
    pca_results.npc=nPCs;
    pca_results.p=p;
else
    pca_results.npc=params_pca.npc;
end
disp(['Number of PCs: ',num2str(pca_results.npc)])

% do the pca
disp('Performing PCA...')
[coeff,score]=fast_pca(X',pca_results.npc);

pc_gene_name(length(hvgix),pca_results.npc)="";
pc_gene_ix(length(hvgix),pca_results.npc)=0;
for i=1:length(coeff(1,:))
    [cs,ix]=sort(coeff(:,i),'descend');
    pc_gene_name(:,i)=gene_name(hvgix(ix));
    pc_gene_id(:,i)=gene_id(hvgix(ix));
    pc_gene_ix(:,i)=hvgix(ix);
end

pca_results.coeff=coeff;
pca_results.coords=score;
pca_results.pc_gene_name=pc_gene_name;
pca_results.pc_gene_id=pc_gene_id;
pca_results.pc_gene_ix=pc_gene_ix;

pcatime=toc;

        
if doPlot
    
    if ~isempty(params_tsne) && numel(params_tsne.initY)==2
        pcix=params_tsne.initY;
    else
        pcix=[1,2];
    end
    
    signix=sign(pcix);
    pcix=abs(pcix);
    S=score(:,pcix);
    if signix(1)<0
        S(:,1)=-S(:,1);
    end
    if signix(2)<0
        S(:,2)=-S(:,2);
    end
    
    figure(2);clf
    plotScatter(S,'group',group,colors,2);
    %     axis square
    drawnow
end


% 3. tSNE on first npc PCs
tic
if ~isempty(params_tsne)
    rng(params_tsne.rngSeed)
    
    if isempty(params_tsne.initY)
        initY=[];
        
    elseif isnumeric(params_tsne.initY)&&numel(params_tsne.initY)==2
        %initY=[pc1,pc2] indicates index of PCs to use as initial points
        
        pcix=params_tsne.initY;
        signix=sign(pcix);
        pcix=abs(pcix);
        initY=score(:,pcix);
        if signix(1)<0
            initY(:,1)=-initY(:,1);
        end
        if signix(2)<0
            initY(:,2)=-initY(:,2);
        end
    
%         initY=score(:,params_tsne.initY);
%         initY=initY+ 5*(rand(size(initY))-0.5); %random jitter?
        
    elseif size(params_tsne.initY,1)==size(score,1)
        initY=params_tsne.initY;
        
    end
    
    disp('Performing tSNE...')
%     res=tsne(score,'perplexity',params_tsne.perplexity,'initialY',initY,...
%         'verbose',2,'NumPrint',100,'options',statset('maxIter',params_tsne.maxiter));
    res=tsne(score,'perplexity',params_tsne.perplexity,'learnRate',params_tsne.learnrate,'exaggeration',params_tsne.exagg,...
        'initialY',initY, 'verbose',2,'NumPrint',100,'options',statset('maxIter',params_tsne.maxiter));
    
    if doPlot
        figure(3);clf
        plotScatter(res,'group',group,colors,3);
        % axis square
        drawnow
    end
    
    tsne_results.coords=res;
    
else
    tsne_results=[];
end
tsnetime=toc;

% 4. Clustering:  could do HC (OK, but optimal leaforder is slow!), kmeans (fast!), densitypeaks, ...
tic
if exist('params_clust','var')&&~isempty(params_clust)
    switch params_clust.method
        case 'kmeans'
            disp('Performing K-means clustering...')
            clusterID=kmeans(pca_results.coords,params_clust.K);
            clust_results.K=params_clust.K;
            
        case 'linkage'
            disp('Performing hierarchical clustering...')
            Z = linkage(pca_results.coords,params_clust.linkage,params_clust.metric);
            clusterID = cluster(Z,'Maxclust',params_clust.K);
            clust_results.Z=Z; %the tree
            clust_results.K=params_clust.K;
            
        case 'modularity'
            disp('Performing graph-based clustering...')
            %Seurat uses KNN=10 and pruning=1/15 as defaults...
            D=knnGraph(pca_results.coords,params_clust.KNN,params_clust.metric);
            COMTY = cluster_jl_cpp(D,1,1,0,0);
            clusterID=COMTY.COM{end};
            clust_results.K=length(COMTY.SIZE{end});
    end
    
    %relabel to be in decreasing size
    nT=arrayfun(@(x)nnz(clusterID==x),1:clust_results.K);
    [nT,ix]=sort(nT,'descend');
    IDold=clusterID;
    for i=1:clust_results.K
        clusterID(IDold==ix(i))=i;
    end
    clust_results.clusterID=clusterID;
 
    if doPlot && ~isempty(tsne_results)
        figure(4);clf
        plotScatter(res,'group',clusterID,colors,4);
        % axis square
        drawnow
    end
    
else
    clust_results=[];
end
clusttime=toc;
% disp('')
disp(['HVG time: ',num2str(hvgtime),', PCA time: ',num2str(pcatime),', tSNE time: ',num2str(tsnetime),', Cluster time: ',num2str(clusttime)])

result.hvg=hvg_results;
result.pca=pca_results;
result.tsne=tsne_results;
result.clust=clust_results;


% n=1:dr.npc;
% %weird way to prepend pc to 1:npc for labels
% pcnames=cellfun(@(x)num2str(x),num2cell(n(:)),'uniformoutput',false);
% pcnames=strcat('pc',pcnames);
% pctable=array2table(score,'variablenames',pcnames);

% dr.tsne.loss=loss;
% n=1:2;
% tsnenames=cellfun(@(x)num2str(x),num2cell(n(:)),'uniformoutput',false);
% tsnenames=strcat('tsne',tsnenames);
% tsnetable=array2table(res,'variablenames',tsnenames);
% dr.data=[pctable,tsnetable];