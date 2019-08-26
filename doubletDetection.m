function [doublets,rep_doublets,rep_counts,tSNE, group]=doubletDetection(rawcounts, params,method)
% inspired by Shor's Doublet detection (python3):
% https://github.com/JonathanShor/DoubletDetection

%my mod: use KNN neighborhoods instead of graph-clustering.
% two main pars: K and Boostrate.
% K - size of kNN nhood. small affects

% several other small mods too
% slow part: resampling parents

%TODO:
% - is FDR really needed?  Also, the cutoff should be fixed given nKnn, nSynth, nCells.
% - options for newLibSize mode (max now, also possible: mean)
% - replace parents true/false (true now)
% - Project synthetic cells onto existing PCs of parent's gene space?


% if ~isempty(params.rngSeed) && isscalar(params.rngSeed)
%     s=RandStream('dsfmt19937','Seed',params.rngSeed);
% else
%     s=RandStream('dsfmt19937');
% end

use_true_cell_HVG=false; %if true, don't recompute HVG on synth-augmented normcounts. saves a little time, doesn't seem to impact negatively
resample_HVG_only=false; %if not, use all genes - pre-feb 15 set this true [probably incorrect]
use_parents_npc=true; %alternative would be to determine new nPC using perm test- very slow.
doMultCompareCorrect=1;

if ~exist('method','var')
    method='knn';
end

if params.nReps==0
    doublets=false(1,nCells);
    rep_counts=0;
    return
end

normcounts=normalizeCounts(rawcounts);

%use true cell HVGs 
if use_true_cell_HVG
    
    %get the HVGs from raw data only.
%     if isfield(params.hvg,'ix') 
%         hvgix=params.hvg.ix; %this can use wrong genes -> genes with zero counts?? pca barfs
%     else
    hvgix=findVariableGenes(normcounts,params.hvg);
%     end
% 
    %subset the counts matrices down to only hvgs to accelerate
    if resample_HVG_only
        rawcounts=rawcounts(hvgix,:);
    end
    
    normcounts=normcounts(hvgix,:);
else
    
end

logcounts=log10(normcounts+1);

[nGenes,nCells]=size(logcounts);

tSNE=[];
group=[];
rep_doublets=false(params.nReps,nCells);

% augment the population with synthetic doublets
% each doublet created from two randomly selected parents => correct proportions
% sum counts - downsample to new lib size.
nSynth=floor(params.boostrate*nCells);

% random #s involved: could repeat several times, and find concensus doublets.
tenth=round(params.nReps/10);

% rawcounts_aug=rawcounts;
% rawcounts_aug(nGenes,nCells+nSynth)=0;
% normcounts_aug=normcounts;
% normcounts_aug(nGenes,nCells+nSynth)=0;
logcounts_aug=logcounts;
logcounts_aug(nGenes,nCells+nSynth)=0;
for it=1:params.nReps
    
%     tic
    rawcounts_synth=createSyntheticDoublets(rawcounts,nSynth);
%     t_gensynth=toc
    
%     tic
%     rawcounts_aug(:,nCells+1:end)=rawcounts_synth;
    normcounts_synth=normalizeCounts(rawcounts_synth);
    if use_true_cell_HVG
        if resample_HVG_only
            logcounts_aug(:,nCells+1:end)=log10(normcounts_synth+1);
        else
            logcounts_aug(:,nCells+1:end)=log10(normcounts_synth(hvgix,:)+1);
        end
    else
        hvgix=findVariableGenes([normcounts,normcounts_synth],params.hvg);
        logcounts_aug=[logcounts(hvgix,:),log10(normcounts_synth(hvgix,:)+1)];
    end

    X=logcounts_aug;
% t_normloghvg=toc
    
%     normcounts_synth=normalizeCounts(rawcounts_synth);
%     
%     if use_true_cell_HVG
%         logcounts_true=logcounts;
%     else
%         hvgix=findVariableGenes([normcounts,normcounts_synth],params.hvg);
%         logcounts_true=logcounts(hvgix,:);
%     end
%     normcounts_synth=normcounts_synth(hvgix,:);
%    
%     logcounts_synth=log10(normcounts_synth+1); %only the new synthDoublets
%     
%     X=[logcounts_true,logcounts_synth];
    
% tic
    %scaling
    X=(X-mean(X,2))./std(X,[],2); %zscore gene distributions
    X(X>params.pca.maxScaled)=params.pca.maxScaled;
    X(X<-params.pca.maxScaled)=-params.pca.maxScaled; 
    
    %compute truncated PCA as space within which to measure distances
    if use_parents_npc
        npc=params.pca.npc;
        [~,score]=fast_pca(X',npc);
    else
        p=findSignificantPCs(X',params.pca.permute_reps,0.05);
        npc=find(~(p==1/params.pca.permute_reps),1,'first')-1
        [~,score]=fast_pca(X',npc);
    end
% t_pca=toc  
    
    
    this_rep_doublets=false(1,nCells);
    switch method
        
        case 'graph_clust'
            
            %use graph-based clustering to group like cells together, compute
            %fraction of each cluster that is synthetic.
            
            % cluster the augmented dataset, Louvain method (granularity set by K-num k-nearest neighbors)
            D=knnGraph(score,params.nKnn,'jaccard');
            COMTY = cluster_jl_cpp(D,1,1,0,0);
            commID=COMTY.COM{end};
            commCount=COMTY.SIZE{end};
            nClust=length(commCount)
            
            %real cells in clusters with high synthetic doublet fraction are labeled as doublets
            synth(nCells+1:end)=true;
            cellscores=zeros(1,nCells); %score for each true cell
            trueCellID=commID(1:nCells);
            synthCount=zeros(nClust,1);
            synthFrac=zeros(nClust,1);
            for j=1:nClust
                synthCount(j)=nnz(synth(commID==j));
                synthFrac(j)=synthCount(j)/commCount(j);
                cellscores(trueCellID==j)=synthFrac(j);
            end
            
            %largest jump heuristic
            %             synthFracSorted=sort(synthFrac,'descend');
            %             [~,cix]=max(abs(diff(synthFracSorted)));
            % cutoff=synthFracSorted(cix);
            % doublets(cellscores>=cutoff)=true;
            % doubletClusters=synthFrac>=cutoff;
            
            % Use significance test for deciding which are doublet clusters:
            %p-values using hypergeometric test: what is probability of observed number of synth or greater drawn from community of
            %given size, given that there are total M cells of which N are synthetic doublets.
            p=hygecdf(synthCount(:),nCells+nSynth,nSynth,commCount(:),'upper');
            %correct for number of clusters
            [h,crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,params.FDR);
            cutoff=min(synthFrac(h));
            this_rep_doublets(cellscores>=cutoff)=true;
            % doubletClusters=synthFrac>=cutoff;
            
            % isequal(doubletClusters(:),h(:))
            
        case 'knn'
            
            %get the kNN of each cell, compute fraction of synth doublets in
            %this neighborhood.
            
            %synthCount=number of true cell's neighbors that are synth
            %synthFrac=fraction of true cell's neighbors that are synth
            
            %reps - keep track of mean synthCount per cell across nReps?
%             tic
            nnix=knnsearch(score,score,'K',params.nKnn+1); nnix(:,1)=[];
            synthCount=sum(nnix(1:nCells,:)>nCells,2); %all synthetic cells have index >nCells
            %             synthFrac=synthCount/params.nKnn;
            p=hygecdf(synthCount(:),nCells+nSynth,nSynth,params.nKnn,'upper');
            
            if doMultCompareCorrect
                [h,crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,params.FDR);
%                 this_rep_doublets=h;
                this_rep_doublets=adj_p<params.FDR;
            else
                this_rep_doublets=p<params.pThr;    
            end
%             t_knn=toc
            
            %             cutoff=min(synthFrac(h));
            %             this_rep_doublets(synthFrac>=cutoff)=true;
    end
    
    rep_doublets(it,:)=this_rep_doublets;
    
    
    group=zeros(1,size(score,1));
    group(1:nCells)=this_rep_doublets; %trudoubs called
    group(nCells+1:end)=2; %synth
    group=categorical(group,0:2,{'tc','td','sd'});

    if isfield(params,'doplot')&&params.doplot
        
%         initY=score(:,params.tsne.initY);

        initY=params.tsne.initY;
        initY=[initY;min(initY)+(max(initY)-min(initY)).*rand(nSynth,2)];
        
%         tSNE=tsne(score);
        tSNE=tsne(score,'perplexity',params.tsne.perplexity,'learnRate',params.tsne.learnrate,'exaggeration',params.tsne.exagg,...
            'initialY',initY, 'verbose',2,'NumPrint',100,'options',statset('maxIter',params.tsne.maxiter));
        plotScatter(tSNE,'group',group,[],it+41);
        drawnow
    end
    
    if tenth>0 && mod(it,tenth)==0
        fprintf('.')
    elseif tenth==0
        fprintf('.')
    end
end
fprintf('\n')

rep_counts=sum(rep_doublets,2);

%voting methods for calling doublets accounting for all rep results
switch params.votemethod
    case 'any'
        % returns all doublets encountered in any rep
        doublets=any(rep_doublets,1);
    case 'all'
        % returns only consensus doublets that were called in all reps
        doublets=all(rep_doublets,1);
    case 'majority'
        % returns doublets that were called in at least half of reps
        doublets=sum(rep_doublets,1)>=0.5*params.nReps;
    case 'p10'
        % returns doublets that were called in more than 10% of reps
        doublets=sum(rep_doublets,1)>0.1*params.nReps;
    case 'p90'
        % returns doublets that were called in more than 90% of reps
        doublets=sum(rep_doublets,1)>0.9*params.nReps;
    otherwise
        warning('invalid option for aggregating doublets found across repetitions')
        doublets=all(rep_doublets,1);
end

end

function [rawcounts_synth,parentIDs,max_parent,synth_libsize]=createSyntheticDoublets(rawcounts,nSynth)

[nGenes,nCells]=size(rawcounts);
rawcounts_synth=zeros(nGenes,nSynth);
parentIDs=zeros(2,nSynth);
max_parent=zeros(1,nSynth);
synth_libsize=zeros(1,nSynth);
for j=1:nSynth
    parentIDs(:,j)=randi(nCells,1,2);
    parents=rawcounts(:,parentIDs(:,j)); %sample parents with replacement
    
    parentLibSize=sum(parents,1);
    [synth_libsize(j),ix]=max(parentLibSize); %assume doublets have more reads?
    %         newLibSize=round(mean(parentLibSize)); %options: max/min/mean/etc?
    max_parent(j)=parentIDs(ix,j);
    
    newSynth=sum(parents,2);
    edges=[0;cumsum(newSynth)];
    perm=randperm(sum(parentLibSize),synth_libsize(j)); %sample gene counts without replacement
    rawcounts_synth(:,j)=histcounts(perm,edges)';
%     perm=randperm(sum(parentLibSize));
%     synthDoublets(:,j)=histcounts(perm(1:newLibSize),edges)';
end


    
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