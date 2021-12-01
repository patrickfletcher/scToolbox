function [ax,hs,cb]=plotDotPlot(varNames,groupNames,sizeData,colorData,figOrAxis,varGroup,sortmethod,params)
arguments
    varNames
    groupNames
    sizeData
    colorData
    figOrAxis = []
    varGroup = []
    sortmethod = 'none'
    params = []
end
%dot plot: area = %>0, color = mean expression value, for a list of genes and
%groups of cells.
%
% [ax,hs,cb]=plotDotPlot(varNames,groupNames,sizeData,colorData,figOrAxis,varGroup,sortmethod,sp_params,do_var_norm)

%TODO: pass this in as varNames
% varLabels=strcat(repmat({'\it '},size(varNames)),varNames);

%TODO: varGroup contains label info, add as text centered on each block

%TODO: pass all the following in via sp_params
% - alt: args block

%TODO: wrapper that operates on ncounts + genes? then could give just names
%list, cellgroup, options. e.g. log=T/F, mean_only_expressed, ...

mrkr='o'; %squares don't scale right

def_params.gap=0.1;
def_params.marg_h=[0.1,0.1];
def_params.marg_w=[0.1,0.1];
def_params.minArea=1;
def_params.maxArea=144;
def_params.min_prct = 5; 
def_params.prct_leg=3; 
def_params.prct_leg_height=0.05;
def_params.cb_prct_gap=0.01;
def_params.cb_gap=0.02;
def_params.cb_width=0.02;
def_params.cblabel='mean log_{10} expression';
def_params.do_var_norm=false;

if isempty(params)
    params=def_params;
end
parfields=fieldnames(def_params);
for i=1:length(parfields)
    if ~isfield(params,parfields{i})
        params.(parfields{i})=def_params.(parfields{i});
    end
end

if isempty(varGroup)
    varGroup=ones(size(varNames));
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
varCats=categories(varGroup);
nPerGroup=countcats(varGroup);

if params.do_var_norm
    colorData=colorData./max(colorData,[],2);
end

switch sortmethod
    case 'alpha'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            varsub=varNames(thiscat);
            sizesub=sizeData(thiscat,:);
            colorsub=colorData(thiscat,:);
            [~,ixs]=natsort(cellstr(varsub));
            varsub=varsub(ixs);
            sizesub=sizesub(ixs,:);
            colorsub=colorsub(ixs,:);
            varNames(thiscat)=varsub;
            sizeData(thiscat,:)=sizesub; 
            colorData(thiscat,:)=colorsub;                  
        end

    case 'size'
        [sizeData,ixs]=groupCountMatrix(sizeData',varGroup,'optim');
        sizeData=sizeData';
        colorData=colorData(ixs,:);
        varNames=varNames(ixs);

    case 'maxsize'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            varsub=varNames(thiscat);
            sizesub=sizeData(thiscat,:);
            colorsub=colorData(thiscat,:);
            [~,ixs]=sort(max(sizesub,[],2),'descend');
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub(ixs,:); 
            sizeData(thiscat,:)=sizesub(ixs,:);                 
        end

    case 'color'
        [colorData,ixs]=groupCountMatrix(colorData',varGroup,'optim');
        colorData=colorData';
        sizeData=sizeData(ixs,:);
        varNames=varNames(ixs);

    case 'maxcolor'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            varsub=varNames(thiscat);
            colorsub=colorData(thiscat,:);
            sizesub=sizeData(thiscat,:);
            [~,ixs]=sort(max(colorsub,[],2),'descend'); 
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub(ixs,:); 
            sizeData(thiscat,:)=sizesub(ixs,:);         
        end
    case 'first'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>params.min_prct); 
            varsub=varNames(thiscat);
            colorsub=colorData(thiscat,:);
            sizesub=sizeData(thiscat,:);
            ixs=unique(r,'stable');
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub(ixs,:); 
            sizeData(thiscat,:)=sizesub(ixs,:);   
        end
    case 'firstsize'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>params.min_prct); 
            varsub=varNames(thiscat);
            colorsub=colorData(thiscat,:);
            sizesub=sizeData(thiscat,:);
            [ixs,ia,ic]=unique(r,'stable');
            varsub=varsub(ixs);
            colorsub=colorsub(ixs,:);
            sizesub=sizesub(ixs,:);

            [sizesub,ixs]=groupCountMatrix(sizesub',c(ia),'optim');
            sizesub=sizesub';
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub(ixs,:); 
            sizeData(thiscat,:)=sizesub;
        end
    case 'firstcolor'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>params.min_prct); 
            varsub=varNames(thiscat);
            colorsub=colorData(thiscat,:);
            sizesub=sizeData(thiscat,:);
            [ixs,ia,ic]=unique(r,'stable');
            varsub=varsub(ixs);
            colorsub=colorsub(ixs,:);
            sizesub=sizesub(ixs,:);

            [colorsub,ixs]=groupCountMatrix(colorsub',c(ia),'optim');
            colorsub=colorsub';
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub; 
            sizeData(thiscat,:)=sizesub(ixs,:);
        end
    otherwise
        %no sorting
end



%%%% rescale frac to [minArea,maxArea]
minArea=params.minArea; 
maxArea=params.maxArea;
% maxArea=(ax(1).Position(3)/length(groupNames))^2;
if minArea==0, minArea=eps; end
if params.min_prct==0, params.min_prct=eps; end

minSizeData=min(sizeData(:));
maxSizeData=max(sizeData(:));
sizeData(sizeData<params.min_prct)=nan;
areaData=rescale(sizeData,minArea,maxArea,'InputMin',minSizeData,'InputMax',maxSizeData);

%prct_leg = # sizes to show >=2 [min%, max% + even breaks between]
if isscalar(params.prct_leg)
    n_prct_leg = params.prct_leg;
    leg_prct=linspace(params.min_prct,maxSizeData,n_prct_leg);
else
    leg_prct=params.prct_leg;
    leg_prct(leg_prct>maxSizeData)=maxSizeData;
    leg_prct=unique(leg_prct);
end
prct_leg_sizes=rescale(leg_prct,minArea,maxArea,'InputMin',minSizeData,'InputMax',maxSizeData);
% prct_leg_sizes=prct_leg_sizes.^2;

prct_leg_nums=round(leg_prct);
prct_leg_labels=strcat(num2str(prct_leg_nums(:)),{'%'});

%%%% set up the figure and axes %%%%%%%%%%%%%%%%
% alt here: use graphics containers...
makeAxes=true;
if isempty(figOrAxis)
    figH=figure();
elseif isa(figOrAxis,'matlab.ui.Figure')
    figH=figure(figOrAxis);clf
elseif isa(figOrAxis,'matlab.graphics.axis.Axes')
    ax=figOrAxis;
    axes(ax);
    makeAxes=false;
end

pos=[params.marg_w(1),params.marg_h(1),1-params.marg_w(1)-params.marg_w(2),1-params.marg_h(1)-params.marg_h(2)];
if makeAxes
%     ax=tight_subplot(1,1,1,sp_params.gap,sp_params.marg_h,sp_params.marg_w);
    ax=axes('OuterPosition',[0,0,1,1],'Position',pos);
end

ax.PositionConstraint='outerposition';
% ax.Position
% ax.PositionConstraint
% ax.OuterPosition

%%%% plot %%%%%%%%%%%%%%%%
[X,Y]=meshgrid(1:length(groupNames),1:length(varNames));
hs=scatter(ax,X(:),Y(:),areaData(:),colorData(:),mrkr,'filled');
hs.MarkerEdgeColor=0.5*[1,1,1];
hs.LineWidth=0.5;

cmap=cbrewer('seq','Reds',64);
% cmap=cbrewer('seq','OrRd',64);
cmap=[1,1,1;cmap];
colormap(ax, cmap)

xlim([0.5,length(groupNames)+0.5])
ylim([0.5,length(varNames)+0.5])
box on

%boost nG separations
if length(nPerGroup)>1
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
cRight=axRight-axPos(3)/2;

cb=colorbar('orientation','horizontal');
cb.Position(2)=axTop+params.cb_gap;
cb.Position(4)=params.cb_width;
cb.Position(3)=axPos(3)/2;
cb.Position(1)=ax.Position(1); %doing this last makes it work...
rp=1;
% cb.Ticks=[ceil(10^rp*cb.Limits(1)),...
%           floor(10^rp*cb.Limits(2))]/10^rp; %round to 1 decimal point
%           round(10^rp*(cb.Limits(1)+diff(cb.Limits)/2)),...
% cb.Limits=cb.Ticks([1,3]);
cb.TickLength=0.025;

cb.Label.String=params.cblabel;
% cb.Label.HorizontalAlignment='right';
cb.Label.Units='normalized';
% cb.Label.Position(1)=0.25;

%legend for %
ax(2)=axes(gcf,'Position',...
    [cRight+params.cb_prct_gap, axTop+params.cb_gap,...
     axRight-cRight-params.cb_prct_gap, params.prct_leg_height]);
 
xvals=ones(size(prct_leg_labels));
yvals=1:length(prct_leg_labels);
hs(2)=scatter(xvals,yvals,prct_leg_sizes,'w',mrkr,'filled');
hs(2).ZData=-prct_leg_sizes;
hs(2).MarkerEdgeColor=0.5*[1,1,1];
hs(2).LineWidth=0.5;
ht=text(ax(2),xvals+0.5,yvals,prct_leg_labels,'HorizontalAlignment','left','FontSize',cb.Label.FontSize);
xlim([0.5,2.75])
ylim([0.5,length(prct_leg_labels)+0.5])
axis off

dcm_obj = datacursormode(figH);
set(dcm_obj,'UpdateFcn',{@customDataTips,groupNames,varNames,sizeData(:),colorData(:)})

end

function txt = customDataTips(~,event_obj,groupNames,varNames,Sz,Col)
% Customizes text of data tips:
%    dcm_obj = datacursormode(fig);
%    set(dcm_obj,'UpdateFcn',{@addGeneDatatips,genes.name})

pos = get(event_obj,'Position');
idx = get(event_obj,'DataIndex');
txt = {[varNames{pos(2)},', ',groupNames{pos(1)}],...
       ['Expr: ',num2str(Col(idx))],...
       ['Prct: ',num2str(Sz(idx))]};
end