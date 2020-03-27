% function [ax,c,lCG,htypelabs]=plotHeatMap(X,geneNames,ident,groupColors,subclusterOpt,figID,nPerGroup,doCellTypeBars,doTypeLabels,nsubsample,PCscores)
function [ax,c,lCG,htypelabs,subsamp]=plotHeatMap(X,varNames,obsGroup,obsGroupColors,figID,varGroup,nsubsample,doTypeLabels,subclusterOpt,sp_params,PCscores)

% heatmap: cells vs genes. genes labeled. cells(genes) optionally grouped by categories, specified in categorical arrays

%bar_width=zero => no bars, empty obsNames => no text labels

%TODO: clean up interface & options
%todo: manage groups of zero size consistently

%%%% input parsing
% p=inputParser();

if ~exist('doTypeLabels','var')
    doTypeLabels=true;
end
% doCellTypeLegend=false;

%unpack inputs
% cellGroupAxLabel='cell type';
%could pass cell array of names in groupBy - use first for grouping, plot remaining as secondary marker(s)

cmap=flipud(cbrewer('div','RdBu',15));
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

if ~exist('subclusterOpt')
    subclusterOpt='none';
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


if ~iscategorical(obsGroup)
    obsGroup=categorical(obsGroup);
end
groupNames=categories(obsGroup);

if ~exist('varGroup','var') || isempty(varGroup)
    varGroup=ones(size(varNames));
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
nPerGroup=countcats(varGroup);


%subsample if requested
if ~exist('nsubsample','var')
    nsubsample=0;
end
subsamp=[];
if nsubsample>0
    for i=1:length(groupNames)
        ix=find(obsGroup==groupNames(i));
        if length(ix)>nsubsample
            thissub=ix(randi(length(ix),nsubsample,1));
        else
            thissub=ix;
        end
        subsamp=[subsamp,thissub];
    end
    X=X(:,subsamp);
    obsGroup=obsGroup(subsamp);
end

groupCounts=countcats(obsGroup);

if ~exist('obsGroupColors','var')||size(obsGroupColors,1)<length(groupNames)
    warning('number of colors does not match number of groups')
    obsGroupColors=cbrewer('qual','Set1',length(groupNames));
    %alt: just cycle back through colors?
end

% HMgeneAxlabel='Genes';
HMgeneAxlabel='Marker genes';
cbLabel='log_{10} fold threshold';
geneLabels=strcat(repmat({'\it '},size(varNames)),varNames);



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
axHM=axes('Position',[sp_params.marg_w(1),sp_params.marg_h(1),...
    1-2*sum(sp_params.marg_w),1-2*sum(sp_params.marg_h)-bar_gap-bar_leg_width]);
% axHM=axes('Position',[sp_params.marg_w(1),sp_params.marg_h(1),...
%     1-sum(sp_params.marg_w)-cb_gap-cb_width,1-sum(sp_params.marg_h)-bar_gap-bar_leg_width]);

%%%% build the grouped count data matrix
if exist('PCscores','var')
[X,cellPerm]=groupCountMatrix(X,obsGroup,subclusterOpt,PCscores);
else
[X,cellPerm]=groupCountMatrix(X,obsGroup,subclusterOpt);
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

imagesc(axHM,X);
colormap(axHM,cmap);

%boost nG separations
yGridValues=cumsum(nPerGroup)+0.5;
line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues(:)','color',0.25*[1,1,1])

set(axHM,'xtick',cvals(2,keepticks(1:end-1))+0.5,'xticklabel',[])
axHM.XGrid='on';
axHM.GridColor=0.3*[1,1,1];
axHM.GridAlpha=.7;

set(axHM,'ytick',1:length(geneLabels),'yticklabel',geneLabels)
% ylabel(axHM,HMgeneAxlabel)

% axHM.Units='points';
axPos=axHM.Position;


%%%% colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cbPos=axPos;
cbPos(1)=axPos(1)+axPos(3)+cb_gap;
cbPos(3)=cb_width;
cbPos(2)=axPos(2)+axPos(4)/6;
cbPos(4)=2*axPos(4)/3;

% c=colorbar(axHM);
% c.Position=cbPos;
c=colorbar(axHM,'Position',cbPos);

% cy=ylabel(c,cbLabel);
% cy.Position(1)=cy.Position(1)*0.66;

%             c.Limits=[0,max(colors(:,i))];
if any(X(:)<0)
    c.Ticks=[ceil(10*c.Limits(1)),0,floor(10*c.Limits(2))]/10; %round to 1 decimal point
    CLIM=caxis(axHM);
    cExt=max(abs(CLIM));
    axHM.CLim=[-cExt,cExt]; %symmetric color axis so middle color of divergent colormaps is at zero
    ylim(c,CLIM)
else
    c.Ticks=[0,floor(10*c.Limits(2))]/10; %round to 1 decimal point
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
        lCG(i).Color=obsGroupColors(i,:);
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
    htypelabs=text(middles,barTextOffset*ones(size(middles)),groupNames(groupCounts>0),...
        'horizontalalignment','center','verticalalignment','bottom','tag','typelab');
end

%%%% secondary axis for geneGroupMarkers (and/or gene dendrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax(3)=axes('position', POS2);

%%%% final cleanup
axHM.YColor=[0,0,0];
linkaxes([axHM,axBL],'x');
axis(axHM,'tight');
axis(axBL,'off')
ax=[axHM,axBL];