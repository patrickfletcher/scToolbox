function [ax,hs,cb]=plotDotPlot(varNames,groupNames,sizeData,colorData,figOrAxis,varGroup,sortmethod,sp_params,do_var_norm)
%dot plot: area = %>0, color = mean expression value, for a list of genes and
%groups of cells.

%TODO> BUG:  ax(2) prct legend sizes are different than main fig (different data axes?)

%%%% input parsing
% p=inputParser();

%TODO: pass this in as varNames
% varLabels=strcat(repmat({'\it '},size(varNames)),varNames);

%TODO: varGroup contains label info, add as text centered on each block

%TODO: pass all the following in via sp_params (or default constructor?)
mrkr='s';
%TODO: finer control of labels for % etc

if ~exist('sp_params','var')||isempty(sp_params)
    sp_params.gap=0.1;
    sp_params.marg_h=0.1;
    sp_params.marg_w=0.1;
    sp_params.minArea=0;
    sp_params.maxArea=200;
    sp_params.prct_leg=25:25:100; %should these be generated from data?
    sp_params.cb_prct_gap=0.01;
    sp_params.cb_gap=0.02;
    sp_params.cb_width=0.02;
    sp_params.cblabel='mean log_{10} expression';
end

if ~exist('do_var_norm','var')
    do_var_norm=false;
end

minArea=sp_params.minArea; 
maxArea=sp_params.maxArea;
% maxArea=(ax(1).Position(3)/length(groupNames))^2;
if minArea==0, minArea=eps; end
prct_leg_areas=minArea+(maxArea-minArea)*sp_params.prct_leg/100;
prct_leg_areas=prct_leg_areas.^2;
prct_leg_labels=strcat(num2str(sp_params.prct_leg(:)),{'%'});

% cmap=cbrewer('qual','Set1',9);
% cmap=cbrewer('seq','Greens',15);
cmap=cbrewer('seq','Reds',15);
% cmap=cmap(5:end,:);

if ~exist('varGroup','var') || isempty(varGroup)
    varGroup=ones(size(varNames));
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
nPerGroup=countcats(varGroup);


if exist('sortmethod','var') && ~isempty(sortmethod)
    switch sortmethod
        case 'alpha'
            [varNames,ixs]=sort(varNames);
            colorData=colorData(ixs,:);
            sizeData=sizeData(ixs,:);
%             varLabels=varLabels(ixs);
        case 'size'
            [sizeData,ixs]=groupCountMatrix(sizeData',varGroup,'optim');
            sizeData=sizeData';
            colorData=colorData(ixs,:);
            varNames=varNames(ixs);
%             varLabels=varLabels(ixs);
        case 'color'
            [colorData,ixs]=groupCountMatrix(colorData',varGroup,'optim');
            colorData=colorData';
            sizeData=sizeData(ixs,:);
            varNames=varNames(ixs);
%             varLabels=varLabels(ixs);
    end

end


[X,Y]=meshgrid(1:length(groupNames),1:length(varNames));

%%%% set up the figure and axes
mustMakeAxes=true;
if exist('figOrAxis','var')
    if isa(figOrAxis,'matlab.ui.Figure')
        figH=figure(figOrAxis);clf
    elseif isa(figOrAxis,'matlab.graphics.axis.Axes')
        ax=figOrAxis;
        axes(ax);
        mustMakeAxes=false;
    end
else
    figH=figure();
end

if mustMakeAxes
    ax=tight_subplot(1,1,1,sp_params.gap,sp_params.marg_h,sp_params.marg_w);
end

colormap(cmap)

%rescale frac to [minArea,maxArea]

%can't have markers of size zero
sizeData=(sizeData-min(sizeData(:)))./(max(sizeData(:))-min(sizeData(:))); %force sizedata into [0,1]
sizeData=minArea+(maxArea-minArea)*sizeData; %area represents prct: 
sizeData=sizeData.^2;

%quantize the size data into bins of 5%
% ds=(maxArea-minArea)/20;
% sizeData=ds*floor(sizeData/ds);
% sizeData(sizeData==0)=eps;

%meanExpr row-wise unitize?
if do_var_norm
    colorData=colorData./max(colorData,[],2);
end

hs=scatter(ax,X(:),Y(:),sizeData(:),colorData(:),mrkr,'filled');
hs.MarkerEdgeColor=0.5*[1,1,1];

xlim([0.5,length(groupNames)+0.5])
ylim([0.5,length(varNames)+0.5])
box on

%boost nG separations
if exist('nPerGroup','var') && ~isempty(nPerGroup)
    nPerGroup=nPerGroup(:)';
    yGridValues=cumsum(nPerGroup(1:end-1))+0.5;
    line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues,'color',0.*[1,1,1],'tag','groupdiv','linewidth',0.5)
end

ax.YDir='reverse';
set(ax,'XTick',1:length(groupNames),'XTickLabel',groupNames)
set(ax,'YTick',1:length(varNames),'YTickLabel',varNames)
axPos=ax.Position;
axTop=axPos(2)+axPos(4);
axRight=axPos(1)+axPos(3);

cb=colorbar('orientation','horizontal');
rp=2;
cb.Ticks=[ceil(10^rp*cb.Limits(1)),floor(10^rp*diff(cb.Limits)/2),floor(10^rp*cb.Limits(2))]/10^rp; %round to 1 decimal point
cb.TickLength=0.025;
ylabel(cb,sp_params.cblabel)

cb.Position(1)=axPos(1);
cb.Position(2)=axTop+sp_params.cb_gap;
cb.Position(3)=axPos(3)/2;
cb.Position(4)=sp_params.cb_width;
cRight=axRight-0.2;

%legend for %
ax(2)=axes(gcf,'Position',...
    [cRight+sp_params.cb_prct_gap, axTop+sp_params.cb_gap,...
     axRight-cRight-sp_params.cb_prct_gap, 1-axTop-sp_params.cb_gap-sp_params.cb_width]);
 
xvals=ones(size(prct_leg_labels));
yvals=1:length(prct_leg_labels);
hs(2)=scatter(xvals,yvals,prct_leg_areas,'w',mrkr,'filled');
hs(2).ZData=prct_leg_areas;
hs(2).MarkerEdgeColor=0.5*[1,1,1];
ht=text(ax(2),xvals+0.75,yvals,prct_leg_labels,'HorizontalAlignment','left','FontSize',cb.Label.FontSize);
xlim([0.5,2.75])
ylim([0.5,length(prct_leg_labels)+0.5])
axis off