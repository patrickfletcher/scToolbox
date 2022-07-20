function result=doubletDetection(rawcounts, normcounts, logcounts, optin, params, hvgpars, pcapars,umappars, options)
arguments
    rawcounts
    normcounts
    logcounts
    
    %optional input data
    optin.hvg=[] 
    optin.pca=[]
    optin.genes=[]
    optin.clusterID=[]

    params.nK=round(logspace(log10(3),log10(sqrt(nnz(cg.subset))),10));
    params.pN=0.15:0.05:0.3
    params.rngSeed=42
    params.pExp=0.075 %how many doublets to remove - reduce to account for homotypic

    params.do_new_pca=1
    params.nn_metric='correlation'

    params.do_new_hvg=0

    params.select_parent_method='random' %implement cluster based
    params.libsize_method='add'
    params.libsize_factor=1 %if using libsize_method='foldmax'

    params.call_method="topk"
    params.do_p_correction=true
    params.pthr=0.05

    hvgpars.nHVG=5000

    pcapars.npc=50
    pcapars.verbose=false
    
    %umap metric and knn will be slaved to params.nn_metric
    umappars.n_neighbors = 30
    umappars.min_dist = 0.3
    umappars.n_epochs = 0
    umappars.n_components=2

    options.top_k = 1
    options.store_all=false
    options.figID=[]
end
% inspired by Shor's Doublet detection (python3):
% https://github.com/JonathanShor/DoubletDetection
% also - scrublet, doubletFinder, scDblFinder... knn-methods.

%TODO: two use cases: 1) param sweep, 2) single run (with optional plotting). split out core work function?
%TODO: option to reuse sweep results but call different expected #

%TODO: implement do_new_hvg

%TODO: interface - linked options (hvgix vs normcounts)

%TODO: -use clustering to select parents from different clusters (homotypic
%doublets don't help). While loop: discard homotypic, continue until nSynth

rng(params.rngSeed,'simdTwister') %for speedup?

[nGenes,nCells]=size(rawcounts);
nExp=round(params.pExp*nCells);

NK=params.nK;
PN=params.pN;

%initialize results
result = params;

%get the HVGs from raw data only
hvgargs=namedargs2cell(hvgpars);
hvg=optin.hvg;
if isempty(hvg) || ~isfield(hvg,'ix')
    hvg=findVariableGenes(normcounts, optin.genes, hvgargs{:}); 
end
% if ~params.do_new_hvg
% nGenes=length(hvg.ix);
% logcounts=logcounts(hvg.ix,:);
% end

%get the true data PCA space (coeffs for projection of new data)
pca=optin.pca;
if isempty(pca)
    pcaargs=namedargs2cell(pcapars);
    pca = doPCA(logcounts(hvg.ix,:), optin.genes, pcaargs{:});
end
pca_score = pca.coords;

% %expand logcounts for synth cells
% % nSynth=floor(params.boost_rate*nCells);
% nSynth=floor(nCells/(1-params.boost_rate) - nCells); %doubletFinder method
% logcounts_aug=logcounts;
% logcounts_aug(nGenes,nCells+nSynth)=0; 
% pca_score(nCells+nSynth, pca.npc)=0; 

%interface: two cases. 1) param sweep, 2) final run supporting plot
% - put the main work in functions (loop body)
% - looping over one PN and one PK naturally works. use "squeeze" to
% flatten after if needed.

% different data structure is needed for param sweep...
% - K by BR grid of synthScores (3D array)
% - for each k, get BC
% - BCmvn=mean(BC)/var(BC), mean&var across BR.
% - find peak BCmvn (save this curve, for verfication)
% -- also save the hygecdf vals...
doublet_scores=zeros(length(NK), nCells, length(PN));
p_hyge=ones(length(NK), nCells, length(PN));
doublets = false(length(NK), nCells, length(PN));
% doublets_hyge = false(length(PK), nCells, length(PN));
BC=zeros(length(PN),length(NK));

% tick=round(params.nReps/10); %progress meter
% t0=tic;
% rep_score=zeros(params.nReps,nCells);
% rep_doublet=false(params.nReps,nCells);
% rep_p=ones(params.nReps,nCells);
% rep_adj_p=ones(params.nReps,nCells);
% for it=1:params.nReps
for i=1:length(PN)
    disp("pN = " + string(PN(i)))
    tic

    %expand logcounts for synth cells
    nSynth=floor(PN(i)*nCells);
%     nSynth=floor(nCells/(1-PN(i)) - nCells); %doubletFinder method
    if params.do_new_hvg
        normcounts_aug=normcounts;
        normcounts_aug(nGenes,nCells+nSynth)=0; 
    end
    logcounts_aug=logcounts;
    logcounts_aug(nGenes,nCells+nSynth)=0; 
    pca_score = pca.coords;
    pca_score(nCells+nSynth, pca.npc)=0; 

    for j=1:length(NK)

        rawcounts_synth=createSyntheticDoublets(rawcounts,nSynth,params, optin.clusterID);
    
        %normalize synthetic doublets across all genes (not just hvgs)
        synth_ncounts=normalizeCounts(rawcounts_synth);

        if params.do_new_hvg
            normcounts_aug(:,nCells+1:end)=synth_ncounts;
            hvg=findVariableGenes(normcounts_aug, optin.genes, hvgargs{:}); 
        end

        synth_lcounts=log10(synth_ncounts(hvg.ix,:)+1);
        logcounts_aug(hvg.ix,nCells+1:end)=synth_lcounts;
        
        if params.do_new_pca
            pcaargs=namedargs2cell(pcapars);
            newpca = doPCA(logcounts_aug(hvg.ix,:), optin.genes, pcaargs{:});
            pca_score = newpca.coords;
        else
            %project synthetic doublets into PC space
            X=synth_lcounts;
            X=pca.scale_fun(X, pca.scale_args);
            X(X>pca.maxScaled)=pca.maxScaled;
            X(X<-pca.maxScaled)=-pca.maxScaled;
            pca_score(nCells+1:end,:)=X'*pca.coeff;
        end
                
        %get the kNN of each cell, compute fraction of synth doublets in
        %this neighborhood.
%         K = round(PK(j)*nCells);
        K = NK(j);
        [nnix, nndists]=knnsearch(pca_score, pca_score,K=K+1,Distance=params.nn_metric); 
        nnix(:,1)=[]; nndists(:,1)=[]; 

        %synthCount=number of true cell's neighbors that are synth
        synthCount=sum(nnix(1:nCells,:)>nCells,2); %all synthetic cells have index >nCells
        synthScore=synthCount/K;
        p=hygecdf(synthCount(:),nCells+nSynth, nSynth, K,'upper');
%         [~,~, ~, adj_p]=fdr_bh(p, params.p_thr);

        doublet_scores(j,:,i)=synthScore;
        p_hyge(j,:,i)=p;

        %call doublets in two ways: hyge and top n (expected)
        switch params.call_method
            case "topk"
                [~, top_ix]=maxk(synthScore,nExp);
                doublets(j,top_ix,i)=true;
            case "pval"
                if params.do_p_correction
                    [~,~,~,p]=fdr_bh(p);
                end
                doublets(j,p<params.pthr,i)=true;
        end
        
        fprintf('.')
    end
    BC(i,:)=bimodality_coefficient(doublet_scores(:,:,i),2); %BC computed across cells

    fprintf('(%0.1fs)',toc)
    fprintf('\n')
end

meanBC=mean(BC,1);
varBC=var(BC,0,1);
if length(PN)==1
    varBC=ones(size(meanBC));
end
BCmvn = meanBC./varBC;
[topBCmvn, topBCKix]=maxk(BCmvn, options.top_k);

%use the top K selection to choose the top N
for i=1:options.top_k
    [b, ix]=max(BC(:,topBCKix(i)));
    topBC(i)=b;
    topBCNix(i)=ix;
end
%alt: could do "voting" on the results from best NK. eg. keep any cell
%called as doublet across those reps (makes sense if same PN used for
%several reps)
% usedoublets=sum(doublets(topBCKix,:,:),3)/length(PN)>=votethr;

usedoublets=squeeze(doublets(topBCKix,:,topBCNix));
usescores=squeeze(doublet_scores(topBCKix,:,topBCNix));
usep=squeeze(p_hyge(topBCKix,:,topBCNix));

result.bestNK=NK(topBCKix);
result.bestPN=PN(topBCNix);
result.doublets=usedoublets;
result.doublet_scores=usescores;
result.p_hyge=usep;

result.BC=BC;
result.BCmvn=BCmvn;
result.topBCmvn=topBCmvn;
result.topBC=topBC;

%for detailed investigation:
if options.store_all
result.all_scores=doublet_scores;
result.all_p_hyge=p_hyge;
result.all_doublets=doublets;
end

% rep_counts=sum(rep_doublet,2);
% 
% %voting methods for calling doublets accounting for all rep results
% if params.vote_thr==0
%     doublet=sum(rep_doublet,1)>params.vote_thr*params.nReps;
% else
%     doublet=sum(rep_doublet,1)>=params.vote_thr*params.nReps;
% end

% result.doublet=doublet;
% result.rep_score=rep_score;
% result.rep_p=rep_p;
% result.rep_adj_p=rep_adj_p;
% result.rep_doublet=rep_doublet;
% result.rep_counts=rep_counts;
% result.nSynth=nSynth;

%do a umap for the final rep for demo purposes
umap=[];
group=[];
groupcols=[0.8,0.8,0.8; 0,0.8,0; 0,0,0];
if ~isempty(options.figID)

    group=zeros(1,size(pca_score,1));
    group(1:nCells)=usedoublets;
    group(nCells+1:end)=2; %synth
    group=categorical(group,[0,2,1],{'tc','sd','td'});

    %could just do a new sim at the best K/N combo...

    knn.idx=nnix; knn.dists=nndists;
    umappars.metric = params.nn_metric;
    umapargs=namedargs2cell(umappars);
    umap=doUMAP(pca_score, knn, umapargs{:});

    figure(options.figID); clf
    scatter_grp(umap.coords, group, gcols=groupcols, fig=options.figID);
    axis tight equal
    drawnow

    result.umap=umap;
    result.group=group;
    result.groupcols=groupcols;
end

end

function rawcounts_synth=createSyntheticDoublets(rawcounts, nSynth, params, clustID)

[nGenes,nCells]=size(rawcounts);
rawcounts_synth=zeros(nGenes,nSynth);
for j=1:nSynth
    switch params.select_parent_method
        case 'random'
            parents=rawcounts(:,randi(nCells,1,2)); %sample parents with replacement
        case 'cluster'
            %ensure parents are from different clusters (force heterotypic parents)
            %select based on cluster proportions
%             tbl=tabulate(clusterID);
%             prop_c=tbl.Count./length(clusterID);
            pix1=randi(nCells); %first parent select randomly
            p1clust=clustID(pix1);
            remaining=find(clustID~=p1clust);
            pix2=remaining(randi(length(remaining)));
            parents=rawcounts(:,[pix1,pix2]); 
    end
    
    parentLibSize=sum(parents,1); %options: max/min/mean/etc? 
    %options for new libsize: max/mean/other? assume doublets have more reads?
    doSubSample=true;
    newSynth=sum(parents,2);
    switch params.libsize_method
        case 'add'
            doSubSample=false;
        case 'simplemean'
            newSynth=newSynth/2;
            doSubSample=false;
        case 'simplemax'
            newSynth=max(parents,[],2);
            doSubSample=false;
        case 'mean'
            synth_libsize=round(mean(parentLibSize));
        case 'min'
            synth_libsize=min(parentLibSize); 
        case 'max'
            synth_libsize=max(parentLibSize);  
        case 'foldmax'
            %assume doublets have more reads?
            synth_libsize=round(max(parentLibSize)*params.libsize_factor);
    end

    if doSubSample
        edges=[0;cumsum(newSynth)];
        perm=randperm(sum(parentLibSize),synth_libsize); %sample gene counts without replacement
        newSynth=histcounts(perm,edges)';
    end
    rawcounts_synth(:,j)=newSynth;
end
    
end