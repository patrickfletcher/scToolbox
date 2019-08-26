function highDispIx=findVariableGenes(ncounts,params,figID,includeNames,gene_ix)
%operates on normalized counts, which maintains mean-variance relations (negative binomial)

%TODO: minimum % expressing instead of params.minExpr? use # of cells?

% if ~exist('params.nBins','var') || isempty(params.nBins)
%     params.nBins=20;
% end
% 
% if ~exist('params.minExpr','var') || isempty(params.minExpr) || params.minExpr==0 
%     params.minExpr=eps;
% end

doPlot=exist('figID','var');


%normalized dispersion:
cellsExpr=sum(ncounts>0,2);
fracExpr=cellsExpr/size(ncounts,2);
meanExpr=mean(ncounts,2,'omitnan');
stdExpr=std(ncounts,[],2,'omitnan');
varExpr=stdExpr.^2;
% dispersion=sqrt(stdExpr);
% dispersion=stdExpr;
% dispersion=varExpr;
% dispersion=varExpr./meanExpr;
% dispersion=(varExpr-meanExpr)./meanExpr;
dispersion=(varExpr-meanExpr)./meanExpr.^2; %negative if mean>var

quantiles=prctile(meanExpr,0:round(100/params.nBins):100);
quantiles(end+1)=max(meanExpr);
quantiles = unique(quantiles);
params.nBins=length(quantiles);

if params.nBins<2
    normdispersion=dispersion;
    
else
    
    Y=discretize(meanExpr,quantiles);
    
    means=zeros(params.nBins,1); %mean
    stds=zeros(params.nBins,1); %standard deviation
    meds=zeros(params.nBins,1); %median
    mads=zeros(params.nBins,1); %median absolute deviation
    % highDispIx=[];
    normdispersion=zeros(params.nBins,1);
    for i=1:params.nBins
        thisIx=find(Y==i);
        thisDisp=dispersion(Y==i);
        
        if nnz(thisIx)>0
            
            switch params.normMethod
                case 'zscore'
                    %Macoscko 2015
                    means(i)=mean(thisDisp);
                    if nnz(thisIx)>1 %to be able to compute variance
                        stds(i)=std(thisDisp);
                    else
                        stds(i)=1; %is this correct? makes the normDisp=0
                    end
                    thisDispZ=(thisDisp-means(i))./stds(i);
                    
                case 'medmad'
                    %Zheng2017 (10X)
                    meds(i)=median(thisDisp);
                    if nnz(thisIx)>1 %to be able to compute variance
                        mads(i)=mad(thisDisp,1); %median absolute deviation
                    else
                        mads(i)=1;
                    end
                    factor=1.49; %to approximate units of stdev (https://www.ibm.com/support/knowledgecenter/en/SS4QC9/com.ibm.solutions.wa_an_overview.2.0.0.doc/modified_z.html)
                    if mads(i)==0
                        factor=1.25; 
                        mads(i)=mad(thisDisp); %mean absolute deviation
                    end
                    thisDispZ=(thisDisp-meds(i))./(mads(i)*factor);
            end
            
            normdispersion(thisIx)=thisDispZ(:);
        end
    end
    
end

%remove genes below params.minExpr from consideration:
normdispersion(meanExpr<params.minExpr)=nan;
normdispersion(cellsExpr<params.minCells)=nan;

%remove negative dispersion genes?
normdispersion(varExpr<meanExpr)=nan;

switch params.selectMethod
    case 'number'
        nHVG=params.nHVG;
        %get the top nHVG genes
        [~,ixs]=sort(normdispersion,'descend','MissingPlacement','last');
        highDispIx=ixs(1:nHVG);
        
    case 'threshold'
        thresh=params.dispThr;
        highDispIx=find(normdispersion>thresh);
        
    case 'allValid'
        highDispIx=find(~isnan(normdispersion));
        
    otherwise
        error('Unknown option')
end


if doPlot
    figure(figID);clf
    
    y=dispersion;
%     y=normdispersion;
%     y=normdispersion+min(normdispersion);

%     y=y+abs(min(y));
%     ycenter=abs(min(y))*[1,1];

    %plot only genes with non-negative dispersion (to avoid plot's warnings in log scale)
%     nonnegative=y>=0;
    nonnegative=true(size(y));
    
    highdisp=false(size(y));
    highdisp(highDispIx)=true;
    
    plot(meanExpr(~highdisp&nonnegative),y(~highdisp&nonnegative),'k.')
    hold on
    plot(meanExpr(highdisp&nonnegative),y(highdisp&nonnegative),'r.');
    
%     if ~isempty(ycenter)
%         plot(xlim,ycenter,'k--')
%     end
    
    if exist('includeNames','var')&&~isempty(includeNames)
        highlightIx=gene_ix;
        line(meanExpr(highlightIx),y(highlightIx),'color','r','marker','o','linestyle','none');
        text(meanExpr(highlightIx),y(highlightIx)+0.03,includeNames,'color',[0.8,0,0]);
    end
    
    
    xlabel('mean expression')
    ylabel('dispersion')
    
    set(gca,'xscale','log','yscale','log')
    
end