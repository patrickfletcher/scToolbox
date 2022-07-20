function result=findVariableGenes(ncounts, genes, block, params, options)
arguments
    ncounts
    genes
    block = []
    params.nBins=25
    params.minExpr=1e-6
    params.minCells=1
    params.center="median"
    params.scale="mad1"
    params.selectMethod='number'
    params.nHVG=3000
    params.dispThr=1.3

    params.combine_val="disp"
    params.combine_method='mean'
    params.equiweight=true

    params.factor_mad0=1.2533;
    params.factor_mad1=1.4826;

    options.all_stats=false
    options.figID=[]
    options.highlightGenes=[]
    options.show_normDisp=0
end
%operates on normalized counts, which maintains mean-variance relations (negative binomial)

%factors to approximate units of stdev (https://www.ibm.com/docs/en/cognos-analytics/11.1.0?topic=terms-modified-z-score)

%TODO: clean up 
% - options.center, options.scale
% - splitapply on bins?

%TODO: add blocking capability
% - combineVar (scran) averages the statistics across batches:
% --> measure normdisp for each batch, then average, then rank.

%plotting
%TODO: tooltips for gene-name
%TODO: switch for plot disp or norm disp

disp('Identifying highly variable genes...')
tic

result = params;

if isempty(block)
    block=ones(1,size(ncounts,2));
end
if ~iscategorical(block)
    block=categorical(block);
end
blocks=categories(block);
n_blocks=length(blocks);

keeps=false(size(ncounts,1),n_blocks);
ncells=zeros(size(ncounts,1),n_blocks);
expr=zeros(size(ncounts,1),n_blocks);
vars=zeros(size(ncounts,1),n_blocks);
disp0=zeros(size(ncounts,1),n_blocks);
dispZ=zeros(size(ncounts,1),n_blocks);
dispZrank=zeros(size(ncounts,1),n_blocks);
for j=1:n_blocks
    
    this_block=block==blocks{j};
    N=ncounts(:,this_block);
    cellsExpr=sum(N>0,2);
    meanExpr=mean(N,2,'omitnan');
    varExpr=var(N,[],2,'omitnan');
    dispersion=(varExpr-meanExpr)./meanExpr.^2; %negative if mean>var
    
    keep=meanExpr>params.minExpr & cellsExpr>params.minCells;
    
    ptiles=prctile(meanExpr(keep),0:100/(params.nBins-1):100);
    ptiles(end+1)=min(meanExpr);
    ptiles(end+1)=max(meanExpr);
    ptiles = sort(unique(ptiles));
    params.nBins=length(ptiles);
    
    Y=discretize(meanExpr,ptiles);
    
    centers=zeros(params.nBins,1);
    scales=zeros(params.nBins,1); 
    normdispersion=zeros(size(dispersion));
    for i=1:params.nBins
%         this_bin=Y==i;
        this_bin= Y==i & keep; %can i exclude non-keep genes this way?
        this_gix=find(this_bin);
        thisDisp=dispersion(this_bin);
        if nnz(this_gix)>0
            switch params.center
                case 'mean'
                    centers(i)=mean(thisDisp);
                case 'median'
                    centers(i)=median(thisDisp);
            end
            switch params.scale
                case 'std'
                    scales(i)=std(thisDisp);
                    factor=1;
                case 'mad0'
                    scales(i)=mad(thisDisp,0);
                    factor=params.factor_mad0;
                    if scales(i)==0 %try falling back to std
                        scales(i)=std(thisDisp);
                        factor=1;
                    end
                case 'mad1'
                    scales(i)=mad(thisDisp,1);
                    factor=params.factor_mad1;
                    if scales(i)==0 %try falling back to meanAD, std
                        scales(i)=mad(thisDisp,0);
                        factor=params.factor_mad0;
                        if scales(i)==0
                            scales(i)=std(thisDisp);
                            factor=1;
                        end
                    end
            end
    
            if nnz(this_gix)==1
                scales(i)=1; %is this correct?
            end
            thisDispZ=(thisDisp-centers(i))./(scales(i)*factor);
            normdispersion(this_gix)=thisDispZ(:);
        end
    end

    %remove genes below minExpr/minCells from consideration:
%     normdispersion(~keep)=nan; %is this the best way?? 
%     normdispersion(~keep)=min(normdispersion(keep));% alt: minimum observed...

    %remove negative dispersion genes?
    % normdispersion(normdispersion<0)=nan;
    % normdispersion(dispersion<0)=nan;

    keeps(:,j)=keep;
    ncells(:,j)=cellsExpr;
    expr(:,j)=meanExpr;
    vars(:,j)=varExpr;
    disp0(:,j)=dispersion;
    dispZ(:,j)=normdispersion;
    [~,ixs]=sort(normdispersion,'descend','MissingPlacement','last');
    dispZrank(:,j)=ixs;
end

keep_all = all(keeps,2);
switch params.combine_val
    case "disp"
        comb_val = dispZ;
        sortdir='descend';
    case "rank"
        comb_val = dispZrank;
        sortdir='ascend';
end

comb_val(~keep_all,:)=nan;
switch params.combine_method
    case 'mean'
        comb_val = mean(comb_val,2); %grand mean across blocks ,'omitnan'
    case 'median'
        comb_val = median(comb_val,2); 
    case 'min'
        comb_val = min(comb_val,[],2); 
    case 'max'
        comb_val = max(comb_val,[],2);
end

[sorted_comb_val,ixs]=sort(comb_val,sortdir,'MissingPlacement',"last");
        
switch params.selectMethod
    case 'number'
        N=min(params.nHVG,nnz(keep_all));
        if N<params.nHVG, warning('Used fewer than requested nHVG'); end
        hvgix=ixs(1:N);
        result=rmfield(result,'dispThr');
        
    case 'threshold'
        hvgix=ixs(sorted_comb_val>params.dispThr);
        result.nHVG=length(hvgix);

    case 'allValid'
        hvgix=ixs(~isnan(sorted_comb_val));
        result.nHVG=length(hvgix);
        result=rmfield(result,'dispThr');

%     case 'minrank'
%         validrank=all(dispZrank<params.nHVG,2);
%         hvgix=ixs(validrank);
%         result.nHVG=length(hvgix);
%         result=rmfield(result,'dispThr');
        
    otherwise
        error('Unknown option')
end

result.nHVG=length(hvgix);
result.ix=hvgix;
result.name=genes.name(hvgix);
% result.gene_id=genes.id(hvgix);

result.comb_val=comb_val;
result.sorted_comb_val=sorted_comb_val;

if options.all_stats
    result.stats.keep=keeps;
    result.stats.ncells=ncells;
    result.stats.mean=expr;
    result.stats.var=vars;
    result.stats.dispersion=disp0;
    result.stats.norm_disp=dispZ;
    result.stats.norm_disp_combined=sorted_comb_val;
    result.stats.rank=dispZrank;
end

disp(mfilename + " time: " + string(toc) +"s")

if ~isempty(options.figID)

rawnormdisp = dispZ;

    figure(options.figID);clf
    
%     y=dispersion; ylab = 'dispersion';
    y=rawnormdisp; ylab = 'norm dispersion';
%     y=normdispersion; ylab = 'norm dispersion';
%     y=normdispersion+min(normdispersion);

    %plot only genes with non-negative dispersion (to avoid plot's warnings in log scale)
    nonnegative=y>=0;
%     nonnegative=true(size(y));
    
    highdisp=false(size(y));
    highdisp(hvgix)=true;
    
    redix=find(highdisp&nonnegative);
    
    if ~isempty(options.highlightGenes)
        highlightIx=[];
        if isnumeric(options.highlightGenes) && isscalar(options.highlightGenes)
            highlightIx = hvgix(1:options.highlightGenes);
            options.highlightGenes = genes.name(highlightIx);
        elseif isstring(options.highlightGenes) || iscellstr(options.highlightGenes)
            highlightIx=getGeneIndices(options.highlightGenes,genes.name);
        end
        ishv=ismember(highlightIx,redix);
        ishvix=highlightIx(ishv);
        isnothv=~ismember(highlightIx,redix);
        isnothvix=highlightIx(isnothv);
    end
    
    
    gap=0.1; margh=[0.15,0.025];margw=[0.1,0.025];
    tight_subplot(1,2,1, gap, margh, margw);
    plot(meanExpr(~highdisp&nonnegative),y(~highdisp&nonnegative),'k.');
    hold on
    plot(meanExpr(highdisp&nonnegative),y(highdisp&nonnegative),'r.');
    
    
    if ~isempty(options.highlightGenes)
        line(meanExpr(ishvix),y(ishvix),'color',[.5,0,0],'marker','.','markersize',10,'linestyle','none');
        text(meanExpr(ishvix),y(ishvix)+0.03,options.highlightGenes(ishv),'color',[.5,0,0]);
        
        line(meanExpr(isnothvix),y(isnothvix),'color',[.5,.5,.5],'marker','.','markersize',10,'linestyle','none');
        text(meanExpr(isnothvix),y(isnothvix)+0.03,options.highlightGenes(isnothv),'color',[.5,.5,.5]);
    end
    
    xlabel('mean expression')
    ylabel(ylab)
    set(gca,'xscale','log','yscale','log')

    
    tight_subplot(1,2,2, gap, margh, margw);
    plot(cellsExpr(~highdisp&nonnegative),y(~highdisp&nonnegative),'k.');
    hold on
    plot(cellsExpr(highdisp&nonnegative),y(highdisp&nonnegative),'r.');
    
    if ~isempty(options.highlightGenes)
        line(cellsExpr(ishvix),y(ishvix),'color',[.5,0,0],'marker','.','markersize',10,'linestyle','none');
        text(cellsExpr(ishvix),y(ishvix)+0.03,options.highlightGenes(ishv),'color',[.5,0,0]);
        
        line(cellsExpr(isnothvix),y(isnothvix),'color',[.5,.5,.5],'marker','.','markersize',10,'linestyle','none');
        text(cellsExpr(isnothvix),y(isnothvix)+0.03,options.highlightGenes(isnothv),'color',[.5,.5,.5]);
    end
    
    xlabel('cells expressing')
%     ylabel(ylab)
    
    set(gca,'xscale','log','yscale','log')
    
end