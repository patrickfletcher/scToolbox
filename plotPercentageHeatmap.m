function [ax, c, ixs]=plotPercentageHeatmap(PRCT,typeLabels,geneNames,varGroup,figOrAxis,sortmethod,sp_params)

%TODO: group indicator/name; gene group as grouping var


if ~exist('sp_params','var')||isempty(sp_params)
    sp_params.gap=0.1;
    sp_params.marg_h=0.1;
    sp_params.marg_w=0.1;
    sp_params.cb_gap=0.02;
    sp_params.cb_width=0.02;
    
    sp_params.c_ticks=0:25:100;
    sp_params.c_binsize=5;
    sp_params.cblabel='% expressing';
end

%TODO: granularity of colormap?
max_prct=max(PRCT(:));
min_prct=min(PRCT(:));
nbins=floor(floor(max_prct-min_prct)/5);
nskip=3;
cmap=cbrewer('seq','Greens',nbins+nskip);
cmap=[1,1,1;cmap(nskip+2:end,:)]; %+2 because nskip+1+white


[nGenes,nTypes]=size(PRCT);

if ~exist('varGroup','var') || isempty(varGroup)
    varGroup=ones(nGenes,1);
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
nPerGroup=countcats(varGroup);

% if exist('nPerGroup','var') && ~isempty(nPerGroup)
%     geneGroup=[];
%     for i=1:length(nPerGroup)
%         geneGroup=[geneGroup;i*ones(nPerGroup(i),1)];
%     end
% else
%     geneGroup=ones(size(geneNames));
% end

if exist('sortmethod','var')
    switch sortmethod
        case 'alpha'
%             [geneNames,ixs]=sort(geneNames);
            [geneNames,ixs]=natsort(cellstr(geneNames));
            PRCT=PRCT(ixs,:);
        case 'nnz'
            [PRCT,ixs]=groupCountMatrix(PRCT',varGroup,'nnz');
            PRCT=PRCT';
            geneNames=geneNames(ixs);
        case 'mean'
            [PRCT,ixs]=groupCountMatrix(PRCT',varGroup,'mean');
            PRCT=PRCT';
            geneNames=geneNames(ixs);
        case 'sum'
            [PRCT,ixs]=groupCountMatrix(PRCT',varGroup,'sum');
            PRCT=PRCT';
            geneNames=geneNames(ixs);
        case 'optim'
            [PRCT,ixs]=groupCountMatrix(PRCT',varGroup,'optim');
            PRCT=PRCT';
            geneNames=geneNames(ixs);
    end
else
    ixs=[];
end


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

imagesc(ax,PRCT)

%minor grid for visual aid
grid_color=0.9*[1,1,1];
xGridValues=1.5:1:nTypes;
yGridValues=1.5:1:nGenes;
line(repmat(xlim()',1,nGenes-1),[1;1]*yGridValues,'color',grid_color)
line([1;1]*xGridValues,repmat(ylim()',1,nTypes-1),'color',grid_color)

%boost nG separations
if exist('nPerGroup','var') && ~isempty(nPerGroup)
    grid_boost_color=0.*[1,1,1];
    nPerGroup=nPerGroup(:)';
    yGridBoost=cumsum(nPerGroup(1:end-1))+0.5;
    line(repmat(xlim()',1,length(yGridBoost)),[1;1]*yGridBoost,'color',grid_boost_color,'tag','groupdiv','linewidth',0.5)
end

ax.XTick=1:nTypes;
ax.XTickLabel=typeLabels;
ax.XTickLabelRotation=90;

geneLabels=strcat(repmat({'\it '},size(geneNames)),geneNames);
ax.YTick=1:nGenes;
ax.YTickLabel=geneLabels;

ax.TickLength=[0,0]; %remove ticks, rely on grid

colormap(ax,cmap)
ax.CLim=[floor(min_prct/5),floor(max_prct/5)]*5;

axPos=ax.Position;

c=colorbar('orientation','horizontal'); 
ylabel(c,sp_params.cblabel)
c.TickLength=0.025;
c.Ticks=sp_params.c_ticks;

c.Position(1)=axPos(1)+0.1*axPos(3);
c.Position(2)=axPos(2)+axPos(4)+sp_params.cb_gap;
c.Position(3)=axPos(3)-0.2*axPos(3);
c.Position(4)=sp_params.cb_width;