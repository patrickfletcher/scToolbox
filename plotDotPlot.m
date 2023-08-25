function hdp=plotDotPlot(varNames,groupNames,sizeData,colorData,figOrAxis,varGroup,options)
arguments
    varNames
    groupNames
    sizeData
    colorData
    figOrAxis = []
    varGroup = []
    options.sortby = 'none'
    options.gap=0.1;
    options.margins=[0.1,0.1,0.1,0.1]
    options.area_range=[1,144];
    options.min_dot_prct = 5; 
    options.prct_leg=3; 
    options.prct_leg_dims=[0.2,0.075];
    options.cb_prct_gap=0.01;
    options.cb_gap=0.01;
    options.cb_dims=[0.25, 0.01];
    options.cblabel='mean log_{10} expression';
    options.cb_digits=2;
    options.do_var_norm=false;
    options.sortgroups='none';
    options.is_diverging=false
    options.cmap=[];
    options.FontSize=9;
    options.xtickangle=90;
end
%dot plot: area = %>0, color = mean expression value, for a list of genes and
%groups of cells.

mrkr='o'; %squares don't scale right

if isempty(varGroup)
    varGroup=ones(size(varNames));
end
if ~iscategorical(varGroup); varGroup=categorical(varGroup); end
varCats=categories(varGroup);
nPerGroup=countcats(varGroup);

if options.do_var_norm
    colorData=colorData./max(colorData,[],2);
end

varNamesIn=varNames;
IXS=[];
switch options.sortby
    case 'alpha'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            varsub=varNames(thiscat);
            sizesub=sizeData(thiscat,:);
            colorsub=colorData(thiscat,:);
            [~,ixs]=natsort(cellstr(varsub));
            varNames(thiscat)=varsub(ixs);
            sizeData(thiscat,:)=sizesub(ixs,:); 
            colorData(thiscat,:)=colorsub(ixs,:); 
%             IXS=[IXS;ixs(:)+length(IXS)];
        end

    case 'size'
        [sizeData,ixs]=groupCountMatrix(sizeData',varGroup,'optim');
        sizeData=sizeData';
        colorData=colorData(ixs,:);
        varNames=varNames(ixs);
%         IXS=[IXS;ixs(:)];

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
%             IXS=[IXS;ixs(:)+length(IXS)];               
        end

    case 'color'
        [colorData,ixs]=groupCountMatrix(colorData',varGroup,'optim');
        colorData=colorData';
        sizeData=sizeData(ixs,:);
        varNames=varNames(ixs);
%         IXS=[IXS;ixs(:)];

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
%             IXS=[IXS;ixs(:)+length(IXS)];         
        end
    case 'first'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>options.min_dot_prct); 
            varsub=varNames(thiscat);
            colorsub=colorData(thiscat,:);
            sizesub=sizeData(thiscat,:);
            ixs=unique(r,'stable');
            varNames(thiscat)=varsub(ixs);
            colorData(thiscat,:)=colorsub(ixs,:); 
            sizeData(thiscat,:)=sizesub(ixs,:);
%             IXS=[IXS;ixs(:)+length(IXS)];   
        end
    case 'firstsize'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>options.min_dot_prct); 
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
%             IXS=[IXS;ixs(:)+length(IXS)];
        end
    case 'firstcolor'
        for i=1:length(varCats)
            thiscat=varGroup==varCats{i};
            [r,c]=find(sizeData(thiscat,:)>options.min_dot_prct); 
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
%             IXS=[IXS;ixs(:)+length(IXS)];
        end
    otherwise
        %no sorting
end

ggrp=ones(1,length(groupNames));
switch options.sortgroups
    case 'size'
        [sizeData,ixs]=groupCountMatrix(sizeData,ggrp,'optim');
        colorData=colorData(:,ixs);
        groupNames=groupNames(ixs);
    case 'color'
        [colorData,ixs]=groupCountMatrix(colorData,ggrp,'optim');
        sizeData=sizeData(:,ixs);
        groupNames=groupNames(ixs);
    otherwise
end

%%%% rescale frac to [minArea,maxArea]
minArea=options.area_range(1); 
maxArea=options.area_range(2);
% maxArea=(ax(1).Position(3)/length(groupNames))^2;
if minArea==0, minArea=eps; end
if options.min_dot_prct==0, options.min_dot_prct=eps; end

minSizeData=min(sizeData(:));
maxSizeData=max(sizeData(:));
sizeData(sizeData<options.min_dot_prct)=nan;
areaData=rescale(sizeData,minArea,maxArea,'InputMin',minSizeData,'InputMax',maxSizeData);

%prct_leg = # sizes to show >=2 [min%, max% + even breaks between]
if isscalar(options.prct_leg)
    n_prct_leg = options.prct_leg;
    leg_prct=linspace(options.min_dot_prct,maxSizeData,n_prct_leg);
else
    leg_prct=options.prct_leg;
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

marg_h=options.margins(1:2);
marg_w=options.margins(3:4);
pos=[marg_w(1),marg_h(1),1-sum(marg_w),1-sum(marg_h)];
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

%handle colormap/diverging
if isempty(options.cmap)
    if options.is_diverging
        cmap=split_cmap();
    else
        cmap=cbrewer('seq','Reds',64); cmap=[1,1,1;cmap];
    end
end
colormap(ax(1),cmap);

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

ax.FontSize=options.FontSize;
xtickangle(ax,options.xtickangle)

cb=colorbar('orientation','horizontal');
rp=options.cb_digits;
cb.Ticks=[ceil(10^rp*cb.Limits(1)),...
          floor(10^rp*cb.Limits(2))]/10^rp; %round to 1 decimal point
%           round(10^rp*(cb.Limits(1)+diff(cb.Limits)/2)),...
% cb.Limits=cb.Ticks([1,3]);
cb.TickLength=0.025;

cb.Label.String=options.cblabel;
% cb.Label.HorizontalAlignment='right';
cb.Label.Units='normalized';
% cb.Label.Position(1)=0.25;

if options.is_diverging
    ax(1).CLim=max(abs(colorData(:)))*[-1,1];
    cb.Limits=[min(colorData(:)),max(colorData(:))];
end

axPos=ax.Position;
axTop=axPos(2)+axPos(4);
axRight=axPos(1)+axPos(3);

cb.Position(1)=axPos(1); 
cb.Position(2)=axTop+options.cb_gap;
cb.Position(3)=options.cb_dims(1);
cb.Position(4)=options.cb_dims(2);

cRight=axPos(1)+options.cb_dims(1);

%legend for %
prct_leg_pos=[1-marg_w(2)-options.prct_leg_dims(1), axTop+options.cb_gap,...
     options.prct_leg_dims(1), options.prct_leg_dims(2)];
% prct_leg_pos=[cRight+options.cb_prct_gap, axTop+options.cb_gap,...
%      axRight-cRight-options.cb_prct_gap, options.prct_leg_height];
ax(2)=axes(gcf,'Position',prct_leg_pos);
 
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

hdp.ax=ax;
hdp.hs=hs;
hdp.cb=cb;
[~,locB]=ismember(varNames,varNamesIn);
hdp.rowperm=locB(:);
hdp.varNames=varNames;
hdp.options=options;

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