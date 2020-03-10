function result = doPCA(tcounts, genes, params, hvgix, figID, pcix)

verbose = false;
result = params;

if exist('hvgix','var')&&~isempty(hvgix)
    X=tcounts(hvgix,:);
else
    X=tcounts;
end

%is rescaling necessary? covariance vs correlation? logtransformed counts are roughly on same scale, in [0,3]
switch params.scale_method
    case 'zscore'
        % Seurat scaling
        X=(X-mean(X,2))./std(X,[],2); %zscore gene distributions
%         if doPlot
%             figure(10);clf
%             histogram(max(X,[],2))
%             %how many cells are >maxScaled? set maxScaled so that this
%             %number is < min expected size of rare subpopulation?
%             histogram(sum(X>params.maxScaled,2))
%             xlabel('# cells Z-score > max')
%             ylabel('# of genes')
%         end
        X(X>params.maxScaled)=params.maxScaled;
        X(X<-params.maxScaled)=-params.maxScaled;  %
        
    case 'unit'
        %unitize
        X=(X-min(X,[],2))./(max(X,[],2)-min(X,[],2));
        
    case 'unit_center'
        %unitize, then subtract mean
        X=(X-min(X,[],2))./(max(X,[],2)-min(X,[],2));
        X=(X-mean(X,2));

    case 'center'
        X=(X-mean(X,2));
        
    case 'none'
end

if params.npc<1
    % Find number of significant PCs
    if verbose, disp('Performing PCA permutation test...'), end
    p=findSignificantPCs(X',params.permute_reps,0.05);
    nPCs=find(~(p==1/params.permute_reps),1,'first')-1;
    result.npc=nPCs;
    result.p=p;
else
    result.npc=params.npc;
end
if verbose, disp(['Number of PCs: ',num2str(result.npc)]), end

% do the pca
if verbose, disp('Computing PCA...'), end
[coeff,score]=fast_pca(X',result.npc);

pc_gene_name(size(X,1),result.npc)="";
pc_gene_ix(size(X,1),result.npc)=0;
for i=1:length(coeff(1,:))
    [cs,ix]=sort(coeff(:,i),'descend');
    if exist('hvgix','var')&&~isempty(hvgix)
        pc_gene_name(:,i)=genes.name(hvgix(ix));
        pc_gene_id(:,i)=genes.id(hvgix(ix));
        pc_gene_ix(:,i)=hvgix(ix);
    else
        pc_gene_name(:,i)=genes.name(ix);
        pc_gene_id(:,i)=genes.id(ix);
        pc_gene_ix(:,i)=ix;
    end
end

result.coeff=coeff;
result.coords=score;
result.pc_gene_name=pc_gene_name;
result.pc_gene_id=pc_gene_id;
result.pc_gene_ix=pc_gene_ix;
        
if exist('figID','var')
    
    if ~exist('pcix','var') || isempty(pcix)
        pcix=[1,2];
    end
    if numel(pcix)~=2
        error("PCA plot: specify two PCs to plot")
    end
    
    %support reversing the direction of pcs
    signix=sign(pcix);
    pcix=abs(pcix);
    S=score(:,pcix);
    S(:,1)=signix(1)*S(:,1);
    S(:,2)=signix(2)*S(:,2);
    
%     group = ones(1,size(X,2));
%     colors=[];
%     figure(figID);clf
%     plotScatter(S,'group',group,colors,figID);
%     plotScatter(S,'group',group,colors,figID);
    
    genelist = pc_gene_name(1,1:2);
    colors=tcounts(pc_gene_ix(1,1:2),:);
    figure(figID);clf
    plotScatter(S,'value',genelist,colors,figID);
    
    %     axis square
    axis tight
    drawnow
end