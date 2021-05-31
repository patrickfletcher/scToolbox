function result=doubletDetection(rawcounts, normcounts, logcounts, params, genes)
% inspired by Shor's Doublet detection (python3):
% https://github.com/JonathanShor/DoubletDetection


%my mod: use KNN neighborhoods instead of graph-clustering.
% two main pars: K and Boostrate.
% K - size of kNN nhood. small affects

% several other small mods too
% slow part: resampling parents

%TODO:
% - Project synthetic cells onto existing PCs of parent's gene space!!!!!
% - is FDR really needed?  Also, the value of synthCount for cutoff should be fixed given n_neighbors, nSynth, nCells.
% - options for newLibSize mode (max now, also possible: mean)
% - replace parents true/false (true now)

rng(params.rngSeed,'simdTwister') %for speedup?
% if ~isempty(params.rngSeed) && isscalar(params.rngSeed)
%     s=RandStream('dsfmt19937','Seed',params.rngSeed);
% else
%     s=RandStream('dsfmt19937');
% end

nCells=size(rawcounts,2);

if ~isnumeric(params.n_neighbors) && params.n_neighbors=="sqrtN"
    params.n_neighbors=floor(sqrt(nCells));
end

%initialize results
result.params = params;

%early exit option: nReps=0 means nothign to do.
if params.nReps==0
    result.doublets=false(1,nCells);
    result.rep_counts=0;
    return
end

%was true, trying false 4/26/21
use_true_cell_HVG=false; %if true, don't recompute HVG on synth-augmented normcounts. saves a little time, doesn't seem to impact negatively
doMultCompareCorrect=true;

if ~isfield(params,'method')
    params.method='knn';
end
if ~isfield(params,'libsize_method')
    params.libsize_method='max';
end

% augment the population with synthetic doublets
% each doublet created from two randomly selected parents => correct
% proportions sum counts - downsample to new lib size.
nSynth=floor(params.boostrate*nCells);

%use true cell HVGs 

hvg=[];
if use_true_cell_HVG
    %get the HVGs from raw data only
    hvg=findVariableGenes(normcounts, genes, params.hvg);
    hvgix=hvg.ix;
  
    %logcounts can be subset now: constant
    logcounts=logcounts(hvgix,:); 
    logcounts_aug=logcounts;
    nGenes=length(hvgix);
    logcounts_aug(nGenes,nCells+nSynth)=0; %expand logcounts for synth cells
    
end

tick=round(params.nReps/10); %progress meter
t0=tic;
rep_doublet=false(params.nReps,nCells);
rep_p=ones(params.nReps,nCells);
rep_adj_p=ones(params.nReps,nCells);
for it=1:params.nReps
    
    rawcounts_synth=createSyntheticDoublets(rawcounts,nSynth,params);

    %normalize synthetic doublets across all genes (not just hvgs)
    normcounts_synth=normalizeCounts(rawcounts_synth);
    
    if use_true_cell_HVG
        logcounts_aug(:,nCells+1:end)=log10(normcounts_synth(hvgix,:)+1);
    else
        hvg=findVariableGenes([normcounts,normcounts_synth], genes, params.hvg);
        logcounts_aug=[logcounts(hvg.ix,:), log10(normcounts_synth(hvg.ix,:)+1)];
    end

    %If always using true_cell_HVG: compute the scale factors for true data
    % (mean, std across cells, per gene). Use those to simply scale the
    % synth data (no recompute).  Then use PCA done on only true data:
    % project scaled_synth onto true data PCs. ==> eliminate redoing PCA
    [~,pca_score]=fast_pca(logcounts_aug', params.pca.npc, params.pca.maxScaled);
    
    this_rep_doublets=false(1,nCells);
    switch params.method
        
        case 'graph_clust'
            
%             %use graph-based clustering to group like cells together, compute
%             %fraction of each cluster that is synthetic.
%             
%             % cluster the augmented dataset, Louvain method (granularity set by K-num k-nearest neighbors)
%             D=knnGraph(score,params.n_neighbors,'jaccard');
%             COMTY = cluster_jl_cpp(D,1,1,0,0);
%             commID=COMTY.COM{end};
%             commCount=COMTY.SIZE{end};
%             nClust=length(commCount)
%             
%             %real cells in clusters with high synthetic doublet fraction are labeled as doublets
%             synth(nCells+1:end)=true;
%             cellscores=zeros(1,nCells); %score for each true cell
%             trueCellID=commID(1:nCells);
%             synthCount=zeros(nClust,1);
%             synthFrac=zeros(nClust,1);
%             for j=1:nClust
%                 synthCount(j)=nnz(synth(commID==j));
%                 synthFrac(j)=synthCount(j)/commCount(j);
%                 cellscores(trueCellID==j)=synthFrac(j);
%             end
%             
%             %largest jump heuristic
%             %             synthFracSorted=sort(synthFrac,'descend');
%             %             [~,cix]=max(abs(diff(synthFracSorted)));
%             % cutoff=synthFracSorted(cix);
%             % doublets(cellscores>=cutoff)=true;
%             % doubletClusters=synthFrac>=cutoff;
%             
%             % Use significance test for deciding which are doublet clusters:
%             %p-values using hypergeometric test: what is probability of observed number of synth or greater drawn from community of
%             %given size, given that there are total M cells of which N are synthetic doublets.
%             p=hygecdf(synthCount(:),nCells+nSynth,nSynth,commCount(:),'upper');
%             %correct for number of clusters
%             [h,crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,params.FDR);
%             cutoff=min(synthFrac(h));
%             this_rep_doublets(cellscores>=cutoff)=true;
%             % doubletClusters=synthFrac>=cutoff;
%             % isequal(doubletClusters(:),h(:))
            
        case 'knn'
            
            %get the kNN of each cell, compute fraction of synth doublets in
            %this neighborhood.
            %synthCount=number of true cell's neighbors that are synth
            
            nnix=knnsearch(pca_score,pca_score,'K',params.n_neighbors+1); nnix(:,1)=[];
            synthCount=sum(nnix(1:nCells,:)>nCells,2); %all synthetic cells have index >nCells
            p=hygecdf(synthCount(:),nCells+nSynth,nSynth,params.n_neighbors,'upper');
            
            [~,~, ~, adj_p]=fdr_bh(p,params.FDR);
            
            if doMultCompareCorrect
                this_rep_doublets=adj_p<params.FDR;
            else
                this_rep_doublets=p<params.pThr;    
            end
    end
    
%     disp("repdoubs: "+num2str(nnz(this_rep_doublets)))
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
if isfield(params,'doplot')&&params.doplot
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
    
end