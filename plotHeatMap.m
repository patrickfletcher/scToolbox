function [ax,hc,lCG,htypelabs,subsamp, cellPerm, varPerm]=plotHeatMap(X,varNames,groups,group_colors,figID,varGroup,nsubsample,doTypeLabels,subclusterOpt,sp_params,PCscores,varClusterOpt)

% heatmap: cells vs genes. genes labeled. cells(genes) optionally grouped by categories, specified in categorical arrays

%bar_width=zero => no bars, empty obsNames => no text labels

%TODO: redo this object oriented (solves many interface problems)
%TODO: clean up interface & options
%todo: manage groups of zero size consistently
%TODO: column vs row vector group sub

%%%% input parsing
% p=inputParser();



% HMgeneAxlabel='Genes';
HMgeneAxlabel='Marker genes';
cbLabel='log_{10} fold threshold';
varLabels=strcat(repmat({'\it '},size(varNames)),varNames);

%unpack inputs
% cellGroupAxLabel='cell type';
%could pass cell array of names in groupBy - use first for grouping, plot remaining as secondary marker(s)

cmap_HM=flipud(cbrewer('div','RdBu',15));
if ~exist('sp_params','var')||isempty(sp_params)
    sp_params.gap=0.1;
    sp_params.marg_h=0.1;
    sp_params.marg_w=0.1;
    sp_params.cb_gap=0.01;
    sp_params.cb_width=0.01;
    sp_params.bar_gap=0.0;
    sp_params.bar_linewidth=4;
    sp_params.bar_leg_width=0.05;
end

if ~exist('doTypeLabels','var')||isempty(doTypeLabels)
    doTypeLabels=true;
end

if ~exist('subclusterOpt')||isempty(subclusterOpt)
    subclusterOpt='none';
end

if ~exist('varClusterOpt')||isempty(varClusterOpt)
    varClusterOpt='none';
end

cb_gap=sp_params.cb_gap;
cb_width=sp_params.cb_width;
bar_linewidth=sp_params.bar_linewidth; 
bar_gap=sp_params.bar_gap;
bar_leg_width=sp_params.bar_leg_width;

doCellTypeBars=true;
if bar_leg_width==0
    doCellTypeBars=false;
    bar_gap=0;
end

if isempty(groups)
    groups=ones(size(X,1),1);
end

if ~iscategorical(groups)
    groups=categorical(groups);
end
groupNames=categories(groups);
groupCounts=countcats(groups);

autoColors=true;
if ~isempty(group_colors)
    if size(group_colors,1)==length(groupNames)
        group_colors=group_colors(groupCounts>0,:);
        autoColors=false;
    end
end

groups=removecats(groups);
groupNames=categories(groups);
groupCounts=countcats(groups);

if autoColors
    group_colors=cbrewer('qual','Set1',max(length(groupNames),3));
end

if ~exist('varGroup','var') || isempty(varGroup)
    varGroup=ones(size(varNames));
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
varGroupCats=categories(varGroup);
nPerGroup=countcats(varGroup);


%subsample if requested
if ~exist('nsubsample','var')
    nsubsample=0;
end
subsamp=[];
if nsubsample>0
    for i=1:length(groupNames)
        ix=find(groups==groupNames(i));
        if length(ix)>nsubsample
%             thissub=ix(randi(length(ix),nsubsample,1));
            thissub=randsample(ix,nsubsample);
        else
            thissub=ix;
        end
        subsamp=[subsamp(:);thissub(:)];
        groupCounts(i)=length(thissub);
    end
    X=X(:,subsamp);
    groups=groups(subsamp);
end

%%%% set up the figure and axes
if exist('figID','var')
    figH=figure(figID);clf
else
    figH=figure();
end

% options: 
% t/f cell mark, if t: t/f dendro
% t/f gene mark, if t: t/f dendro

cellMarkWidth=6; %units: points - axCG position should use points too.

%main axis
% axHM=axes('Position',[sp_params.marg_w(1),sp_params.marg_h(1),...
%     1-2*sum(sp_params.marg_w),1-2*sum(sp_params.marg_h)-bar_gap-bar_leg_width]);
% axHM=axes('Position',[sp_params.marg_w(1),sp_params.marg_h(1),...
%     1-sum(sp_params.marg_w)-cb_gap-cb_width,1-sum(sp_params.marg_h)-bar_gap-bar_leg_width]);

axHM=tight_subplot(1,1,1,0,sp_params.marg_h,sp_params.marg_w);

%sorting on variables too?
if varClusterOpt=="alpha"
    for i=1:length(varGroupCats)
        thisgroup=varGroup==varGroupCats{i};
        thisNames=varNames(thisgroup);
        thisLabels=varLabels(thisgroup);
        thisX=X(thisgroup,:);
        [~,ixs]=sort(varNames(thisgroup));
        varNames(thisgroup)=thisNames(ixs);
        varLabels(thisgroup)=thisLabels(ixs);
        X(thisgroup,:)=thisX(ixs,:);
    end
else
    varPerm=1:length(varNames);
    [Y,varPerm]=groupCountMatrix(X',varGroup,varClusterOpt);
    X=Y';
    varLabels=varLabels(varPerm);
end

%%%% build the grouped count data matrix
if exist('PCscores','var')&&~isempty(PCscores)
[X,cellPerm]=groupCountMatrix(X,groups,subclusterOpt,PCscores);
else
[X,cellPerm]=groupCountMatrix(X,groups,subclusterOpt);
end


%%%% build the cell and gene group marker arrays
% groupby will be contiguous
csum=cumsum(groupCounts); csum=csum(:)';
cvals(1,:)=[0,csum(1:end-1)];
cvals(2,:)=csum;

% non-contiguous group marking (i.e. secondary marker variable, not grouped by): 
% cvals=getContiguousBlocks(groupIdent);

keepticks=any(cvals>0,1);
cvals=cvals(:,groupCounts>0);
keepticks=keepticks(groupCounts>0);

%%%% main axis for heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

him=imagesc(axHM,X);
colormap(axHM,cmap_HM);

%boost nG separations
yGridValues=cumsum(nPerGroup)+0.5;
line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues(:)','color',0.25*[1,1,1])

axHM.XTick=cvals(2,keepticks(1:end-1))+0.5;
axHM.XTickLabel=[];
axHM.XGrid='on';
axHM.GridColor=0.3*[1,1,1];
axHM.GridAlpha=.7;

axHM.YTick=1:length(varLabels);
axHM.YTickLabel=varLabels;

% axHM.Units='points';
axPos=axHM.Position;


%%%% colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cbPos=axPos;
cbPos(1)=axPos(1)+axPos(3)+cb_gap;
cbPos(3)=cb_width;
cbPos(2)=axPos(2)+axPos(4)/2;
cbPos(4)=axPos(4)/2;

hc=colorbar(axHM,'Position',cbPos);

cy=ylabel(hc,cbLabel);
cy.VerticalAlignment='baseline';

%             c.Limits=[0,max(colors(:,i))];
if any(X(:)<0)
    hc.Ticks=[ceil(10*hc.Limits(1)),0,floor(10*hc.Limits(2))]/10; %round to 1 decimal point
    CLIM=caxis(axHM);
    cExt=max(abs(CLIM));
    axHM.CLim=[-cExt,cExt]; %symmetric color axis so middle color of divergent colormaps is at zero
    ylim(hc,CLIM)
else
    hc.Ticks=[0,floor(10*hc.Limits(2))]/10; %round to 1 decimal point
end



%%%% secondary axis for cellGroupMarkers (and/or cell dendrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barLegendPos=axPos;
barLegendPos(2)=axPos(2)+axPos(4)+bar_gap;
barLegendPos(4)=bar_leg_width;
axBL=axes('position', barLegendPos);

% set(axCG,'xtick',[],'color','none','xcolor','none')
% set(axCG,'ytick',cellMarkVal,'yticklabel',cellGroupAxLabel)
% axCG.YAxisLocation='right';
ylim(axBL,[-.5,0.5]) 

barTextOffset=0;
if doCellTypeBars
    barline_x=cvals+0.5;
    barline_y=zeros(size(cvals));
    lCG=line(axBL,barline_x,barline_y);
    for i=1:length(groupNames(groupCounts>0))
        lCG(i).Color=group_colors(i,:);
        lCG(i).LineWidth=bar_linewidth;
    end
    barTextOffset=0.2;
%     if doCellTypeLegend
%         hleg1=legend(lCG,groupNames(groupCounts>0),'Location','westoutside');
%         hleg1.AutoUpdate='off';
%         ax(2).Position=POS2;
%         legend boxoff
%     end
else
    lCG=[];
end


%legend as groupName text at middle of group
if doTypeLabels
    middles=mean(cvals+0.5,1);
    % middles=middles(groupCounts>0);
%     htypelabs=text(middles,barTextOffset*ones(size(middles)),groupNames(groupCounts>0),...
%         'horizontalalignment','center','verticalalignment','bottom','tag','typelab',...
%         'Interpreter','none');
    htypelabs=text(middles,barTextOffset*ones(size(middles)),groupNames(groupCounts>0),...
        'horizontalalignment','center','verticalalignment','bottom','tag','typelab');
else
    htypelabs=[];
end

%%%% secondary axis for geneGroupMarkers (and/or gene dendrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax(3)=axes('position', POS2);

%%%% final cleanup
axHM.YColor=[0,0,0];
linkaxes([axHM,axBL],'x');
axis(axHM,'tight');
axis(axBL,'off')
ax=[axHM,axBL];