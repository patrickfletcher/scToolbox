function result=findVariableGenes(ncounts,genes,options,plotopts)
arguments
    ncounts
    genes
    
    options.normMethod='medmad'
    options.nBins=25
    options.minExpr=eps
    options.minCells=10
    options.selectMethod='number'
    options.nHVG=3000
    options.dispThr=1.3

    plotopts.figID=[]
    plotopts.highlightGenes=[]
    plotopts.show_normDisp=0
end
%operates on normalized counts, which maintains mean-variance relations (negative binomial)

%TODO: add per-batch capability

result = options;

%plotting
%TODO: tooltips for gene-name
%TODO: switch for plot disp or norm disp

%normalized dispersion:
cellsExpr=sum(ncounts>0,2);
% fracExpr=cellsExpr/size(ncounts,2);
meanExpr=mean(ncounts,2,'omitnan');
% stdExpr=std(ncounts,[],2,'omitnan');
% varExpr=stdExpr.^2;
varExpr=var(ncounts,[],2,'omitnan');
% dispersion=sqrt(stdExpr);
% dispersion=stdExpr;
% dispersion=varExpr;
% dispersion=varExpr./meanExpr;
% dispersion=(varExpr-meanExpr)./meanExpr;
dispersion=(varExpr-meanExpr)./meanExpr.^2; %negative if mean>var

clear ncounts

quantiles=prctile(meanExpr,0:100/(options.nBins):100);
quantiles(end+1)=max(meanExpr);
quantiles = unique(quantiles);
options.nBins=length(quantiles);

if options.nBins<2
    normdispersion=dispersion;
    
else
    
    Y=discretize(meanExpr,quantiles);
    
    means=zeros(options.nBins,1); %mean
    stds=zeros(options.nBins,1); %standard deviation
    meds=zeros(options.nBins,1); %median
    mads=zeros(options.nBins,1); %median absolute deviation
    % highDispIx=[];
    normdispersion=zeros(options.nBins,1);
    for i=1:options.nBins
        thisIx=find(Y==i);
        thisDisp=dispersion(Y==i);
        
        if nnz(thisIx)>0
            
            switch options.normMethod
                case 'zscore'
                    %Macoscko 2015
                    means(i)=mean(thisDisp);
                    if nnz(thisIx)>1 %to be able to compute variance
                        stds(i)=std(thisDisp);
                    else
                        stds(i)=1; %is this correct? makes the normDisp=0
                    end
                    thisDispZ=(thisDisp-means(i))./stds(i);
                
                case 'meanmad'
                    %mean + mean absolute deviation
                    meds(i)=mean(thisDisp);
                    if nnz(thisIx)>1 %to be able to compute variance
                        mads(i)=mad(thisDisp,0); %median absolute deviation
                    else
                        mads(i)=1;
                    end
%                     if mads(i)==0
% %                         factor=1.25; 
% %                         mads(i)=mad(thisDisp); %mean absolute deviation
%                         stds(i)=std(thisDisp);
%                     end
                    thisDispZ=(thisDisp-meds(i))./mads(i);
                    
                case 'medmad'
                    %Zheng2017 (10X)
                    meds(i)=median(thisDisp);
                    if nnz(thisIx)>1 %to be able to compute variance
                        mads(i)=mad(thisDisp,1); %median absolute deviation
                    else
                        mads(i)=1;
                    end
                    factor=1;
%                     factor=1.49; %to approximate units of stdev (https://www.ibm.com/support/knowledgecenter/en/SS4QC9/com.ibm.solutions.wa_an_overview.2.0.0.doc/modified_z.html)
                    if mads(i)==0
%                         factor=1.25; 
                        mads(i)=mad(thisDisp); %mean absolute deviation
%                         mads(i)=std(thisDisp);
                    end
                    thisDispZ=(thisDisp-meds(i))./(mads(i)*factor);
            end
            
            normdispersion(thisIx)=thisDispZ(:);
        end
    end
    
end

rawnormdisp = normdispersion;

%remove genes below params.minExpr from consideration:
normdispersion(meanExpr<options.minExpr)=nan;
normdispersion(cellsExpr<options.minCells)=nan;

%remove negative dispersion genes?
% normdispersion(normdispersion<0)=nan;
% normdispersion(dispersion<0)=nan;
% normdispersion(varExpr<meanExpr)=nan;

[sortedNormDisp,ixs]=sort(normdispersion,'descend','MissingPlacement','last');
        
switch options.selectMethod
    case 'number'
        hvgix=ixs(1:options.nHVG);
        result=rmfield(result,'dispThr');
        
    case 'threshold'
        hvgix=ixs(sortedNormDisp>options.dispThr);
        result.nHVG=length(hvgix);

    case 'allValid'
        hvgix=ixs(~isnan(sortedNormDisp));
        result.nHVG=length(hvgix);
        
    otherwise
        error('Unknown option')
end

result.nHVG=length(hvgix);
result.ix=hvgix;
result.name=genes.name(hvgix);
% result.gene_id=genes.id(hvgix);

ishvg=false(size(meanExpr));
ishvg(hvgix)=true;


if ~isempty(plotopts.figID)
    figure(plotopts.figID);clf
    
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
    
    if ~isempty(plotopts.highlightGenes)
        highlightIx=[];
        if isnumeric(plotopts.highlightGenes) && isscalar(plotopts.highlightGenes)
            highlightIx = hvgix(1:plotopts.highlightGenes);
            plotopts.highlightGenes = genes.name(highlightIx);
        elseif isstring(plotopts.highlightGenes) || iscellstr(plotopts.highlightGenes)
            highlightIx=getGeneIndices(plotopts.highlightGenes,genes.name);
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
    
    
    if ~isempty(plotopts.highlightGenes)
        line(meanExpr(ishvix),y(ishvix),'color',[.5,0,0],'marker','.','markersize',10,'linestyle','none');
        text(meanExpr(ishvix),y(ishvix)+0.03,plotopts.highlightGenes(ishv),'color',[.5,0,0]);
        
        line(meanExpr(isnothvix),y(isnothvix),'color',[.5,.5,.5],'marker','.','markersize',10,'linestyle','none');
        text(meanExpr(isnothvix),y(isnothvix)+0.03,plotopts.highlightGenes(isnothv),'color',[.5,.5,.5]);
    end
    
    xlabel('mean expression')
    ylabel(ylab)
    set(gca,'xscale','log','yscale','log')

    
    tight_subplot(1,2,2, gap, margh, margw);
    plot(cellsExpr(~highdisp&nonnegative),y(~highdisp&nonnegative),'k.');
    hold on
    plot(cellsExpr(highdisp&nonnegative),y(highdisp&nonnegative),'r.');
    
    if ~isempty(plotopts.highlightGenes)
        line(cellsExpr(ishvix),y(ishvix),'color',[.5,0,0],'marker','.','markersize',10,'linestyle','none');
        text(cellsExpr(ishvix),y(ishvix)+0.03,plotopts.highlightGenes(ishv),'color',[.5,0,0]);
        
        line(cellsExpr(isnothvix),y(isnothvix),'color',[.5,.5,.5],'marker','.','markersize',10,'linestyle','none');
        text(cellsExpr(isnothvix),y(isnothvix)+0.03,plotopts.highlightGenes(isnothv),'color',[.5,.5,.5]);
    end
    
    xlabel('cells expressing')
%     ylabel(ylab)
    
    set(gca,'xscale','log','yscale','log')
    
end