function [ax, c]=plotExprHeatmap(Expr,typeLabels,geneNames,nPerGroup,figID,sortmethod)

%TODO: random subsamble of K cells

%TODO: group indicator/name; gene group as grouping var

cbLabel='mean normalized counts';
% cbLabel='% max norm counts';

%TODO: granularity of colormap?
cmap=cbrewer('seq','Blues',20);
cmap=[1,1,1;cmap(2:end,:)];


geneLabels=strcat(repmat({'\it '},size(geneNames)),geneNames);
nGenes=length(geneNames);

% if ~iscategorical(geneGroup)
%     geneGroup=categorical(geneGroup);
% end
% nPerGroup=countcats(geneGroup);

% maxExpr=max(Expr,[],2);
% Expr=Expr./maxExpr*100;

if exist('nPerGroup','var') && ~isempty(nPerGroup)
    geneGroup=[];
    for i=1:length(nPerGroup)
        geneGroup=[geneGroup;i*ones(nPerGroup(i),1)];
    end
else
    geneGroup=ones(size(geneNames));
end

if exist('sortmethod','var')
    switch sortmethod
        case 'alpha'
            [geneNames,ixs]=sort(geneNames);
            Expr=Expr(ixs,:);
            geneLabels=geneLabels(ixs);
        case 'optim'
            [Expr,ixs]=groupCountMatrix(Expr',geneGroup,'optim');
            Expr=Expr';
            geneNames=geneNames(ixs);
            geneLabels=geneLabels(ixs);
    end

end



%%%% set up the figure and axes
if exist('figID','var')
    figH=figure(figID);clf
else
    figH=figure();
end

% figH.Position

ax=axes();

colormap(cmap)
imagesc(ax,Expr)

%horizontal lines for visual aid
yGridValues=0.5:1:nGenes;
line(repmat(xlim()',1,nGenes),[1;1]*yGridValues,'color',0.9*[1,1,1])

%boost nG separations
if exist('nPerGroup','var') && ~isempty(nPerGroup)
    nPerGroup=nPerGroup(:)';
    yGridValues=cumsum(nPerGroup(1:end-1))+0.5;
    line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues,'color',0.*[1,1,1],'tag','groupdiv','linewidth',0.5)
end

% for i=1:length(geneLabels)
%     labloc=round(find(Expr(i,:)>maxExpr(i)/10));
%     for j=1:length(labloc)
%         text(labloc(j),i,sprintf('%4.1f',Expr(i,labloc(j))),'Color','w','HorizontalAlignment','center');
%     end
% end

set(ax,'XTick',1:length(typeLabels),'XTickLabel',typeLabels)
set(ax,'YTick',1:length(geneLabels),'YTickLabel',geneLabels)
% ax.Position(1)=ax.Position(1)*0.8;
% ax.Position(2)=ax.Position(2)*0.8;
% ax.Position(3)=ax.Position(3)*1.1;
ax.Position(2)=ax.Position(2)-0.03;
ax.Position(4)=ax.Position(4)-0.03;
axPos=ax.Position;

c=colorbar('orientation','horizontal'); 
ylabel(c,cbLabel)
% c.Ticks=0:25:100;
% c.Ticks=[5,50,100];

% caxis([0,100]);

c.Position(1)=axPos(1)+0.1;
c.Position(2)=axPos(2)+axPos(4)+0.01;
c.Position(3)=axPos(3)-0.2;
c.Position(4)=0.02;