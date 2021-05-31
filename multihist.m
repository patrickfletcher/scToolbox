function [ax,hh,lh]=multihist(X,names,group,showZeros,thresh,colors,figID,panels,sp_params)

%TODO: reorder argument list to make fast use easier (eg. X, thresh, names, group, 
        
%TODO: a stacked bar chart variant

maxBins=30;

doGroupedHistograms=false;
if exist('group','var')&&~isempty(group)
    doGroupedHistograms=true;
    group=group(:)'; %row vector
    if ~iscategorical(group)
        group=categorical(group);
    end
    group=removecats(group);
    groupNames=categories(group);
    typeCounts=countcats(group);

    %don't plot groups with no members
    keep_groups=typeCounts>0;
    nTypes=nnz(keep_groups);
    groupNames=groupNames(keep_groups);
    
    %reorder, smallest to largest counts
    [~,ixs]=sort(typeCounts,'ascend');
    group=reordercats(group,groupNames(ixs));
    groupNames=groupNames(ixs);
    
    if ~exist('colors','var')||size(colors,1)<length(groupNames)
        %make some colors for groups if needed
        colors=cbrewer('qual','Set1',max(length(groupNames),3));
    else
        colors=colors(keep_groups,:);
        colors=colors(ixs,:);
    end
    
    
end

if ~exist('showZeros','var')||isempty(showZeros)
    showZeros=true;
end
if showZeros
    nonzero=true(size(X));
else
    nonzero=X>0;
end

doThresh=false; 
if exist('thresh','var')&&~isempty(thresh)
    doThresh=true;
else
    thresh=0;
end

if isscalar(thresh) && ~isscalar(names)
    thresh=repmat(thresh,1,length(names));
end


if ~exist('figID','var')
    figH=figure();
else
    figH=figure(figID);
end
clf


[nsub,nobs]=size(X);
if nsub>30
    X=X(1:30,:);
    names=names(1:30);
    nsub=30;
end

if ~exist('panels','var')||isempty(panels)
    p=numSubplots(length(names));
    nr=p(1);
    nc=p(2);
else
    nr=panels(1);
    nc=panels(2);
end

if ~exist('sp_params','var')||isempty(sp_params)
%     spmargins=[0.01,0.01];
    sp_params.gap=0.1;
    sp_params.marg_h=0.05;
    sp_params.marg_w=0.05;
end

nBinEdges=min( round(max(10,nobs/10)), maxBins);

for i=1:nsub
    edges=linspace(min(X(i,:)),max(X(i,:)),nBinEdges);
    
    ax(i)=tight_subplot(nr,nc,i,sp_params.gap,sp_params.marg_h,sp_params.marg_w);
    if doGroupedHistograms
        
        %reorder from fewest to most counts of this gene. Need to do
        %different legend for each gene in that case...
        thisGroup = group; thisNames=groupNames; thisCols=colors;
%         meanExpr=zeros(1,nTypes);
%         for j=1:nTypes
%             meanExpr(j)=mean(X(i,thisGroup==thisNames(j) & nonzero(i,:)),2);
%         end
%         [~,ixs]=sort(meanExpr,'descend');
%         thisGroup=reordercats(thisGroup,thisNames(ixs));
%         thisNames=thisNames(ixs);
%         thisCols=thisCols(ixs,:);
        
        hold on
        for j=1:nTypes
            hh(i,j)=histogram(X(i,thisGroup==thisNames(j) & nonzero(i,:)),edges);
            hh(i,j).FaceColor=thisCols(j,:);
        end
    else
        hh(i)=histogram(X(i,nonzero(i,:)),edges,'facecolor',[0.5,0.5,0.5]);

%         histfit(X(i,nonzero(i,:)),nBinEdges-1,'Kernel','Support','positive');
%         histogram(X(i,nonzero(i,:)),edges,'facecolor',[0.5,0.5,0.5],'normalization','probability')
%         [f,xi]=ksdensity(X(i,:));
%         line(xi,f,'color','r')
    end
    
    
    axis tight
    
    if doThresh
        YLIM=ylim();
        line(thresh(i)*[1,1],YLIM,'color','k','linewidth',2)
    end
    
%     xlabel('log_{10} expression')
    title(names(i))
    xlim([min(edges),max(edges)])
end

if doGroupedHistograms
    lh=legend(hh(end,:),groupNames);
    legend boxoff
end

% figH.PaperOrientation='landscape';