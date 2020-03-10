function result=doDimRedAndCluster(scdata,ncounts,tcounts,cells,genes,params, doPlot, group, colors)
% HVGs, PCA, tSNE, UMAP, clustering
% result is a struct with fields hvg, pca, tsne, umap, and clust
% params are copied into the result struct. 

%TODO: Cellranger uses all genes but truncates to fixed # PCs. Here, I use
%fixed number HVG, measure # PCs by perm test... Which is correct?

if ~isempty(scdata) && isstruct(scdata)
    result=scdata; %appends results to existing struct
    ncounts=ncounts(scdata.genesub,scdata.subset);
    tcounts=tcounts(scdata.genesub,scdata.subset);
    genes=genes(scdata.genesub,:);
end

%store the parameters used in the same struct as results
% result.hvg=params.hvg;
% result.pca=params.pca;
% result.tsne=params.tsne;
% result.umap=params.umap;
% result.clust=params.clust;

if ~exist('group','var')||isempty(group)
    group=ones(1,size(ncounts,2));
end
if ~exist('colors','var')||isempty(colors)
    colors=[0.5,0.5,0.5];
end
if ~exist('doPlot','var')||isempty(doPlot)
    doPlot=false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. hvgs
tic
disp('Finding highly variable genes...')
% params.hvg.nHVG=min(round(size(ncounts,2)/2),params.hvg.nHVG);

%can optionally pass extra genes to highlight or force include? store these in
%params?
markers=["POMC","GH1","PRL","LHB","TSHB","SOX2"];
if doPlot
    figure(1);clf
    result.hvg=findVariableGenes(ncounts,genes,params.hvg,1,markers);
    drawnow
else
    result.hvg=findVariableGenes(ncounts,genes,params.hvg);
end
result.hvg.nHVG = length(result.hvg.ix);
result.hvg.nHVG
hvgtime=toc

% for the rest, restrict to HVGs:
X = tcounts(result.hvg.ix,:);
G = genes(result.hvg.ix,:);

% also pick out some top genes to show expression:
meanExpr=mean(ncounts(result.hvg.ix,:),2,'omitnan');
[~,ixm]=sort(meanExpr,'descend','MissingPlacement','last');
valuenames = G.name(ixm(1:4));
colors = X(ixm(1:4),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b. regressout
% X = regressOut([cells.molecPerCell(scdata.subset), cells.genesPerCell(scdata.subset), cells.fracMT(scdata.subset)],X')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. pca
tic
pcplotix = [1,2];
if isfield(params.tsne,'initY') && isnumeric(params.tsne.initY)&&numel(params.tsne.initY)==2
    pcplotix=params.tsne.initY;
end
result.pca = doPCA(X, G, params.pca, [], 2, pcplotix);
disp(['Number of PCs: ',num2str(result.pca.npc)])
pcatime=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. tSNE on first npc PCs
tic
if isfield(params,'tsne')&&~isempty(params.tsne)
    result.tsne = doTSNE(result.pca.coords, params.tsne, 3, valuenames, colors);
else
    result.tsne=[];
end
tsnetime=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. UMAP on first npc PCs
tic
if isfield(params,'umap')&&~isempty(params.umap)
    result.umap = doUMAP(result.pca.coords, params.umap, 4, valuenames, colors);   
else
    result.umap=[];
end
umaptime=toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Clustering:  could do HC (OK, but optimal leaforder is slow!), kmeans (fast!), densitypeaks, ...
tic
if isfield(params,'clust')&&~isempty(params.clust)
 
    result.clust = doClustering(result.pca.coords,params.clust);
    
    if doPlot && ~isempty(result.tsne)
        figure(5);clf
        plotScatter(result.tsne.coords,'group',result.clust.clusterID,colors,5);
        % axis square
        axis tight
        drawnow
    end
    if doPlot && ~isempty(result.umap)
        figure(6);clf
        plotScatter(result.umap.coords,'group',result.clust.clusterID,colors,6);
        % axis square
        axis tight
        drawnow
    end
else
    result.clust=[];
end
clusttime=toc
% disp('')
fmt = '%.2g';
disp(['times - HVG: ',num2str(hvgtime,fmt),', PCA: ',num2str(pcatime,fmt),...
    ', tSNE: ',num2str(tsnetime,fmt),', UMAP: ',num2str(umaptime,fmt),...
    ', Cluster: ',num2str(clusttime,fmt),'s'])


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