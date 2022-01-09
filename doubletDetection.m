function result=doubletDetection(rawcounts, normcounts, logcounts, pca_score, pca_loadings, params)
arguments
    rawcounts
    normcounts
    logcounts
    pca_score
    pca_loadings
    params.rngSeed=42
    params.libsize_method='mean'
    params.use_true_loadings=true
    params.do_p_correction=true
    params.nReps=1
    params.p_thr=0.05
    params.n_neighbors=30
    params.votemethod='any'
    params.doplot=false
    params.hvg=[]
    params.pca=[]
    params.umap=[]
end
% inspired by Shor's Doublet detection (python3):
% https://github.com/JonathanShor/DoubletDetection

%assume PCA loadings + scores of true cells are given. Project synthetic
%doublets into this space.

%my mod: use KNN neighborhoods instead of graph-clustering.
% two main pars: K and Boostrate.
% K - size of kNN nhood. small affects

% several other small mods too
% slow part: resampling parents

%TODO:
% - scanning k for kNN should be much faster than processing different
% boostrates. Store results in a way that allows easy updates?

%TODO:
% - Project synthetic cells onto existing PCs of parent's gene space!!!!!
% - is FDR really needed?  Also, the value of synthCount for cutoff should be fixed given n_neighbors, nSynth, nCells.
% - options for newLibSize mode (max now, also possible: mean)
% - replace parents true/false (true now)

%TODO: -use clustering to select parents from different clusters (homotypic
%doublets don't help). While loop: discard homotypic, continue until nSynth

rng(params.rngSeed,'simdTwister') %for speedup?

nCells=size(rawcounts,2);

%initialize results
result.params = params;

%early exit option: nReps=0 means nothign to do.
if params.nReps==0
    result.doublets=false(1,nCells);
    result.rep_counts=0;
    return
end

% augment the population with synthetic doublets
% each doublet created from two randomly selected parents => correct
% proportions sum counts - downsample to new lib size.
nSynth=floor(params.boostrate*nCells);

% hvg=[];
% if params.use_true_loadings
%     %get the HVGs from raw data only
%     hvg=findVariableGenes(normcounts, genes, params.hvg);
%     hvgix=hvg.ix;
%   
%     %logcounts can be subset now: constant
%     logcounts=logcounts(hvgix,:); 
%     logcounts_aug=logcounts;
%     nGenes=length(hvgix);
%     logcounts_aug(nGenes,nCells+nSynth)=0; %expand logcounts for synth cells
%     
%     %get the true data PCA space (coeffs for projection of new data). Only
%     %possible if use_true_cell_HVG=true
%     [coeff,pca_score, mu, sig]=fast_pca(logcounts', params.pca.npc, params.pca.maxScaled);
%     pca_score(nCells+nSynth,params.pca.npc)=0; 
% end

nGenes=length(hvgix);

logcounts=logcounts(hvgix,:); 
logcounts_aug=logcounts;
logcounts_aug(nGenes,nCells+nSynth)=0; %expand logcounts for synth cells
pca_score(nCells+nSynth,params.pca.npc)=0; 

tick=round(params.nReps/10); %progress meter
t0=tic;
rep_score=zeros(params.nReps,nCells);
rep_doublet=false(params.nReps,nCells);
rep_p=ones(params.nReps,nCells);
rep_adj_p=ones(params.nReps,nCells);
for it=1:params.nReps
    
    rawcounts_synth=createSyntheticDoublets(rawcounts,nSynth,params);

    %normalize synthetic doublets across all genes (not just hvgs)
    normcounts_synth=normalizeCounts(rawcounts_synth);
    
    synth_lcounts=log10(normcounts_synth(hvgix,:)+1);
    logcounts_aug(:,nCells+1:end)=synth_lcounts;
    
    X=synth_lcounts';
    X=(X-mu)./sig;
    X(X>params.pca.maxScaled)=params.pca.maxScaled;
    X(X<-params.pca.maxScaled)=-params.pca.maxScaled;
    pca_score(nCells+1:end,:)=X*coeff;
    
    this_rep_doublets=false(1,nCells);
            
    %get the kNN of each cell, compute fraction of synth doublets in
    %this neighborhood.
    %synthCount=number of true cell's neighbors that are synth
    
    nnix=knnsearch(pca_score,pca_score,'K',params.n_neighbors+1); nnix(:,1)=[];
    synthCount=sum(nnix(1:nCells,:)>nCells,2); %all synthetic cells have index >nCells
    p=hygecdf(synthCount(:),nCells+nSynth,nSynth,params.n_neighbors,'upper');
    
    [~,~, ~, adj_p]=fdr_bh(p,params.p_thr);
    
    if params.do_p_correction
        this_rep_doublets=adj_p<params.p_thr;
    else
        this_rep_doublets=p<params.p_thr;    
    end
    
%     disp("repdoubs: "+num2str(nnz(this_rep_doublets)))
    rep_score(it,:)=synthCount;
    rep_doublet(it,:)=this_rep_doublets; 
    rep_p(it,:)=p;
    rep_adj_p(it,:)=adj_p;  %store p vals so re-thresholding can be done
    
    if tick>0 && mod(it,tick)==0
        fprintf('.%d(%0.1fs)',it,toc)
    elseif tick==0
        fprintf('.')
    end
end
fprintf('\n')

rep_counts=sum(rep_doublet,2);

%voting methods for calling doublets accounting for all rep results
% if params.vote_thr==0
%     doublet=sum(rep_doublet,1)>params.vote_thr*params.nReps;
% else
%     doublet=sum(rep_doublet,1)>=params.vote_thr*params.nReps;
% end

switch params.votemethod
    case 'any'
        % returns all doublets encountered in any rep
        doublet=any(rep_doublet,1);
    case 'all'
        % returns only consensus doublets that were called in all reps
        doublet=all(rep_doublet,1);
    case 'majority'
        % returns doublets that were called in at least half of reps
        doublet=sum(rep_doublet,1)>=0.5*params.nReps;
    case 'p10'
        % returns doublets that were called in more than 10% of reps
        doublet=sum(rep_doublet,1)>0.1*params.nReps;
    case 'p90'
        % returns doublets that were called in more than 90% of reps
        doublet=sum(rep_doublet,1)>0.9*params.nReps;
    otherwise
        warning('invalid option for aggregating doublets found across repetitions')
        doublet=all(rep_doublet,1);
end

%do a umap for the final rep for demo purposes
umap=[];
group=[];
groupcols=[1,1,1; 0,0.8,0; 0,0,0];
if params.doplot
    group=zeros(1,size(pca_score,1));
    group(1:nCells)=doublet;
%     group(1:nCells)=this_rep_doublets; 
    group(nCells+1:end)=2; %synth
    group=categorical(group,[0,2,1],{'tc','sd','td'});

    umap=doUMAP(pca_score, params.umap);
    fh=figure();
    plotScatter(umap.coords, 'group', group, groupcols, fh);
    drawnow
end

result.doublet=doublet;
result.rep_score=rep_score;
result.rep_p=rep_p;
result.rep_adj_p=rep_adj_p;
result.rep_doublet=rep_doublet;
result.rep_counts=rep_counts;
result.nSynth=nSynth;
result.hvg=hvg;
result.umap=umap;
result.group=group;
result.groupcols=groupcols;

end

function rawcounts_synth=createSyntheticDoublets(rawcounts,nSynth,params)

[nGenes,nCells]=size(rawcounts);
rawcounts_synth=zeros(nGenes,nSynth);
for j=1:nSynth
    parents=rawcounts(:,randi(nCells,1,2)); %sample parents with replacement
    
    parentLibSize=sum(parents,1); %options: max/min/mean/etc? 
    %options for new libsize: max/mean/other? assume doublets have more reads?
    switch params.libsize_method
        case 'mean'
            synth_libsize=round(mean(parentLibSize));
        case 'min'
            synth_libsize=min(parentLibSize); 
        case 'max'
            synth_libsize=max(parentLibSize);  
        case 'oversample'
            %assume doublets have more reads?
            synth_libsize=round(max(parentLibSize)*params.libsize_factor);
    end
    
    newSynth=sum(parents,2);
    edges=[0;cumsum(newSynth)];
    perm=randperm(sum(parentLibSize),synth_libsize); %sample gene counts without replacement
    rawcounts_synth(:,j)=histcounts(perm,edges)';
end
    
end




% function [rawcounts_synth,parentIDs,max_parent,synth_libsize]=createSyntheticDoublets(rawcounts,nSynth)
% 
% [nGenes,nCells]=size(rawcounts);
% rawcounts_synth=zeros(nGenes,nSynth);
% parentIDs=zeros(2,nSynth);
% max_parent=zeros(1,nSynth);
% synth_libsize=zeros(1,nSynth);
% for j=1:nSynth
%     parentIDs(:,j)=randi(nCells,1,2);
%     parents=rawcounts(:,parentIDs(:,j)); %sample parents with replacement
%     
%     parentLibSize=sum(parents,1); %options: max/min/mean/etc? 
% %     synth_libsize(j)=round(mean(parentLibSize)); %assume doublets have more reads?
% %     ix=1; %assume doublets have more reads?
%     [synth_libsize(j),ix]=max(parentLibSize); %assume doublets have more reads? 
%     max_parent(j)=parentIDs(ix,j);
%     
%     newSynth=sum(parents,2);
%     edges=[0;cumsum(newSynth)];
%     perm=randperm(sum(parentLibSize),synth_libsize(j)); %sample gene counts without replacement
%     rawcounts_synth(:,j)=histcounts(perm,edges)';
% %     perm=randperm(sum(parentLibSize));
% %     synthDoublets(:,j)=histcounts(perm(1:newLibSize),edges)';
% end


    
    %     %vectorized - different answer, diff rng stream? so far can't vectorize randperm/histcount. is there another way?
    %     parentIDs=randi(nCells,2,nSynth);
    %     parents=cat(3,rawcounts(:,parentIDs(1,:)),rawcounts(:,parentIDs(2,:)));
    %     parentLibSize=squeeze(sum(parents,1))';
    %     newLibSize=max(parentLibSize,[],1);
    %     synthDoublets=zeros(nGenes,nSynth);
    %     for j=1:nSynth
    %         newSynth=sum(parents(:,j,:),3);
    %         edges=[0;cumsum(newSynth)];
    %         perm=randperm(sum(parentLibSize(:,j)));
    %         synthDoublets(:,j)=histcounts(perm(1:newLibSize(j)),edges)';
    %     end
    
    %     rawcounts_synth=[rawcounts(resample_genes,:),synthDoublets];
    %     logcounts_synth=log10(normalizeCounts(rawcounts_synth)+1);