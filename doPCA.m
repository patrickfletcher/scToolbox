function result = doPCA(tcounts, genes, batch, params, options)
arguments
    tcounts
    genes
    batch=[]
    params.hvgix=[]
    params.scale_method='zscore'
    params.center=true
    params.scale=false;
    params.maxScaled=10
    params.npc=0
    params.permute_reps=100
    params.scale_pcs=false
    params.full_coeff=false
    options.verbose=true
    options.figID=[]
    options.pcix=[]
end

%TODO: multibatch PCA - batchelor.  Seems like just a specific choice of
%centering/scaling??  
% - center gene expression using grand mean of batch means
% - scale using weigths
% - project centered back onto PC space computed from scaled..

%TODO: investigate 1) pre-scale/center data, 2) post normalize PCs

%covariance = X'*X/(n-1) <=> X is centered.

% https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca

%TODO: standardized scores = sqrt(n-1)*U
%TODO: clarify loadings vs principal axes (V).  loadings need: *sqrt(S)
%TODO: add interface to numerical options of eigs (tolerance, etc)


result = params;

[nGenes,nCells]=size(tcounts);
if ~isempty(params.hvgix)
    X=tcounts(params.hvgix,:);
else
    X=tcounts;
end

%is rescaling necessary? covariance vs correlation? logtransformed counts are roughly on same scale, in [0,3]

%TODO: create function handle
scale_args = [];
scale_fun = @(x) x;
doscale=false;
switch params.scale_method
    case 'zscore'
        % --> PCA on correlation matrix, not covariance.
        % --> standardized loadings (correlation)

        % Seurat scaling
%         X=(X-mu)./sig; %zscore gene distributions
        
        scale_args.xmean=mean(X,2);
        scale_args.xstd=std(X,[],2);
        scale_fun = @(x, args) (x-args.xmean)./args.xstd;
        doscale=true;
        %notes
        % - largest Zscores occur for small-mean genes.
                
%         if doPlot
%             figure(10);clf
%             histogram(max(X,[],2))
%             %how many cells are >maxScaled? set maxScaled so that this
%             %number is < min expected size of rare subpopulation?
%             histogram(sum(X>params.maxScaled,2))
%             xlabel('# cells Z-score > max')
%             ylabel('# of genes')
%         end
        
    case 'meanmad'
        %use median/mad as more robust measures of center and scale...
        %however, smaller scale values ==> larger standardized. Those will
        %hold more weigh, so clipping even more important...
        MU=mean(X,2);
        MAD=mad(X,0,2); 
        MAD(MU==0 & MAD==0)=1; %not sure if need to have MU=0...???
        
        scale_args.xmean=MU;
        scale_args.xmad=MAD;
        scale_fun = @(x, args) (x-args.xmean)./args.xmad;
        
        
    case 'medmad'
        %use median/mad as more robust measures of center and scale
        MED=median(X,2);
        MAD=mad(X,1,2); 
        MAD(MED==0 & MAD==0)=1; %zero/zero occurs. Setting MAD=1 => X unchanged: X=(X-0)/1
        scale_args.xmed=MED;
        scale_args.xmad=MAD;
        scale_fun = @(x, args) (x-args.xmed)./args.xmad;
        
    case 'unit'
        %unitize
        scale_args.xmin=min(X,[],2);
        scale_args.xmax=min(X,[],2);

        scale_fun = @(x,args) (x-args.xmin)./(args.xmax-args.xmin);
        
    case 'center'
        % --> PCA on covariance matrix
        scale_args = mean(X,2);
        scale_fun = @(x, mu) x-mu;

    case 'none'
        scale_args = [];
        scale_fun = @(x, args) x;
        
end

result.scale_args=scale_args;
result.scale_fun=scale_fun;

X = scale_fun(X, scale_args);

if doscale
X(X>params.maxScaled)=params.maxScaled;
X(X<-params.maxScaled)=-params.maxScaled;
end

if params.npc<1
    % Find number of significant PCs
    if options.verbose, disp('Performing PCA permutation test...'), end
    [npc,p,e,d]=findSignificantPCs(X',params.permute_reps,0.01);
    result.npc=npc;
    result.p=p;
    result.explained=e;
    result.explained_null=d;
    if options.verbose, disp(['Number of PCs: ',num2str(result.npc)]), end
else
    result.npc=params.npc;
end

% do the pca
if options.verbose, disp('Computing PCA...'), end

%useless intermediate function...
% [coeff,score]=fast_pca(X',result.npc);

[U,S,coeff] = svdsecon(X',result.npc);
score =  U*S';

if params.scale_pcs
    score=score./std(score,0,1);
end


% apply the rotation to all genes
% UNTESTED
if params.full_coeff && ~isempty(params.hvgix)
    fullC = zeros(size(tcounts,2),result.npc);
    fullC(params.hvgix,:)=coeff;
    nothvg=setdiff(1:size(tcounts,2), params.hvgix);
    Xleft = tcounts(:, nothvg)';
    Xleft = scale_fun(Xleft, scale_args);
    Xleft = Xleft*U' - mean(Xleft,2)*sum(U,1);
    Xleft = Xleft./diag(S)
end

result.U=U;
result.S=S;
result.coeff=coeff;
result.coords=score;
result.var_exp = diag(S).^2 / (nCells - 1);
result.total_var = sum(var(X,0,2));



% pc_gene_name(size(X,1),result.npc)="";
% % pc_gene_id(size(X,1),result.npc)="";
% pc_gene_ix(size(X,1),result.npc)=0;
% for i=1:length(coeff(1,:))
%     [~,ix]=sort(coeff(:,i),'descend');
%     if ~isempty(params.hvgix)
%         pc_gene_name(:,i)=genes.name(params.hvgix(ix));
% %         pc_gene_id(:,i)=genes.id(params.hvgix(ix));
%         pc_gene_ix(:,i)=hvgix(ix);
%     else
%         pc_gene_name(:,i)=genes.name(ix);
% %         pc_gene_id(:,i)=genes.id(ix);
%         pc_gene_ix(:,i)=ix;
%     end
% end
% result.pc_gene_name=pc_gene_name;
% result.pc_gene_id=pc_gene_id;
% result.pc_gene_ix=pc_gene_ix;
        
if ~isempty(options.figID)

    pc_gene_name(size(X,1),result.npc)="";
    pc_gene_ix(size(X,1),result.npc)=0;
    for i=1:length(coeff(1,:))
        [~,ix]=sort(coeff(:,i),'descend');
        if ~isempty(params.hvgix)
            pc_gene_name(:,i)=genes.name(params.hvgix(ix));
            pc_gene_ix(:,i)=params.hvgix(ix);
        else
            pc_gene_name(:,i)=genes.name(ix);
            pc_gene_ix(:,i)=ix;
        end
    end
    
    if isempty(options.pcix)
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
    
    genelist = [pc_gene_name(1,1:2),pc_gene_name(end,1:2)];
    gix=[pc_gene_ix(1,1:2),pc_gene_ix(end,1:2)];
    vals=tcounts(gix,:);

    figure(options.figID);clf
    hsc=scatter_val(S,vals,fig=options.figID,splitnames=genelist);

    for i=1:length(hsc.ax)
        axis(hsc.ax(i),'equal','tight')
    end
    drawnow
end