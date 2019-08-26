function [ax, c, ixs]=plotPercentageHeatmap(PRCT,typeLabels,geneNames,varGroup,figID,sortmethod,sp_params)

%TODO: group indicator/name; gene group as grouping var


if ~exist('sp_params','var')||isempty(sp_params)
    sp_params.gap=0.1;
    sp_params.marg_h=0.1;
    sp_params.marg_w=0.1;
    sp_params.cb_gap=0.02;
    sp_params.cb_width=0.02;
    
    sp_params.c_ticks=0:25:100;
    sp_params.c_binsize=5;
    sp_params.cbLabel='% expressing';
end

%TODO: granularity of colormap?
nbins=round((max(PRCT(:))-min(PRCT(:)))/5);
cmap=cbrewer('seq','Greens',nbins+3);
cmap=[1,1,1;cmap(5:end,:)];


nGenes=length(geneNames);

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
            [geneNames,ixs]=sort(geneNames);
            PRCT=PRCT(ixs,:);
        case 'optim'
            [PRCT,ixs]=groupCountMatrix(PRCT',varGroup,'optim');
            PRCT=PRCT';
            geneNames=geneNames(ixs);
    end
else
    ixs=[];
end



%%%% set up the figure and axes
if exist('figID','var')
    figH=figure(figID);clf
else
    figH=figure();
end

% figH.Position

ax=tight_subplot(1,1,1,sp_params.gap,sp_params.marg_h,sp_params.marg_w);

colormap(cmap)
imagesc(ax,PRCT)

%horizontal lines for visual aid
yGridValues=0.5:1:nGenes;
line(repmat(xlim()',1,nGenes),[1;1]*yGridValues,'color',0.9*[1,1,1])

%boost nG separations
if exist('nPerGroup','var') && ~isempty(nPerGroup)
    nPerGroup=nPerGroup(:)';
    yGridValues=cumsum(nPerGroup(1:end-1))+0.5;
    line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues,'color',0.*[1,1,1],'tag','groupdiv','linewidth',0.5)
end

set(ax,'XTick',1:length(typeLabels),'XTickLabel',typeLabels)

geneLabels=strcat(repmat({'\it '},size(geneNames)),geneNames);
set(ax,'YTick',1:length(geneLabels),'YTickLabel',geneLabels)
axPos=ax.Position;

c=colorbar('orientation','horizontal'); 
ylabel(c,sp_params.cbLabel)
c.Ticks=0:25:100;
c.TickLength=0.025;

c.Ticks=sp_params.c_ticks;

% caxis([0,100]);
caxis([min(PRCT(:)),max(PRCT(:))]);

c.Position(1)=axPos(1)+0.1*axPos(3);
c.Position(2)=axPos(2)+axPos(4)+sp_params.cb_gap;
c.Position(3)=axPos(3)-0.2*axPos(3);
c.Position(4)=sp_params.cb_width;