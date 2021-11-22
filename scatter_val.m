function hsc=scatter_val(coords,cvals,splitby,opts,scopts)
arguments
    coords
    cvals=[]
    splitby=[]
    opts.title=[] %label to use as title
    opts.fig=[] %specify figure to plot in (if empty, new figure)
    opts.ax=[] %specify axis to plot in (ignored if multipanel)
    opts.margins=[0.1,0.1,0.1,0.1] %[left, bottom, right, top]
    opts.tilegaps=[0.05,0.05] %multi-subplot gaps [horizontal, vertical]
    opts.panels=[] %specific layout of [nrow, ncol] subplots

    opts.cbgap=0.01 %space between ax and cb
    opts.cbdims=[0.3,0.01] %[length, width]
    opts.cblabel=[]
    opts.cbLoc {mustBeMember(opts.cbLoc,["east","west","north","south"])}='east'
    opts.commonCbar=true %true - one cb for all axes same color scale. false - individual cbs
    opts.cmap=[]
    opts.isdiverging=false

    opts.draworder {mustBeMember(opts.draworder,["random","value","flat"])} ='flat'
    scopts.?matlab.graphics.chart.primitive.Scatter
end
%returns a struct containing all handles for post-modification

[nObs,nDim]=size(coords);

% default markersize
% defsz=13-log(nObs);
defsz=5;

%support 3D scatterplots
do3D=false;
scoptargs=namedargs2cell(scopts);
if nDim==2
    scatterfun=@(X,C)scatter(X(:,1),X(:,2),defsz,C,'filled',scoptargs{:});
elseif nDim==3
    scatterfun=@(X,C)scatter3(X(:,1),X(:,2),X(:,3),defsz,C,'filled',scoptargs{:});
    do3D=true;
else
    error('coors are not 2 or 3 dimensional');
end

%default colors
if isempty(cvals)
    cvals=ones(nObs,1);
end

% organize panel splitting
if isempty(splitby)
    splitby=ones(nObs,1);
end
if ~iscategorical(splitby)
    splitby=categorical(splitby);
end
splitby=removecats(splitby);
snames=categories(splitby);
nSplit=length(snames);

%decide whether to make figure and/or axes
doNewAx=false;
if isempty(opts.fig)&&isempty(opts.ax) %make new fig+axes
    fh=figure();
    doNewAx=true;
elseif isempty(opts.fig)&&~isempty(opts.ax) %create correct # ax
    if length(opts.ax)~=nSplit
        error("Cannot fit "+string(nSplit)+" scatterplots into "+string(length(opts.ax))+" axes");
    end
    axes(opts.ax(1))
    fh=gcf;
elseif ~isempty(opts.fig)&&isempty(opts.ax) %make new axes
    fh=figure(opts.fig);
    doNewAx=true;
else %both fig and ax specified: splitby only possible if correct # ax
    if length(opts.ax)~=nSplit
        error("Cannot fit "+string(nSplit)+" scatterplots into "+string(length(opts.ax))+" axes");
    end
    fh=figure(opts.fig);
end

if doNewAx
    nr=floor(sqrt(nSplit));
    nc=ceil(sqrt(nSplit));
    ax=tight_subplot(nr,nc,[],opts.tilegaps,opts.margins([1,3]),opts.margins([2,4]));
end


for i=1:nSplit
    axes(ax(i))

    thisgrp=splitby==snames{i};
    thiscvals=cvals(thisgrp);
    hs(i)=scatterfun(coords(thisgrp,:),thiscvals);
    
    if ~opts.commonCbar
        rectpos=ax(i).Position;
        hcb(i)=makeCB(ax(i),rectpos,opts.cbgap,opts.cbdims,opts.cbLoc);
        if opts.isdiverging
            ax(i).CLim=abs(max(thiscvals))*[-1,1];
%             cbtick=unique([min(thiscvals),0,max(thiscvals)]);
            cbtick=unique(sort([min(hcb(i).Ticks),max(hcb(i).Ticks),0]));
        else
            ax(i).CLim=[min(thiscvals),max(thiscvals)];
%             cbtick=[min(thiscvals),max(thiscvals)];
            cbtick=[min(hcb(i).Ticks),max(hcb(i).Ticks)];
        end
        hcb(i).Limits=[min(thiscvals),max(thiscvals)];
        hcb(i).Ticks=sort(cbtick);
        hcb(i).TickLength=opts.cbdims(2);
    else
        if opts.isdiverging
            ax(i).CLim=abs(max(cvals))*[-1,1];
        else
            ax(i).CLim=[min(cvals),max(cvals)];
        end
    end
    
    axis off
%     axis tight
%     axis equal

    if nSplit>1
        title(snames{i})
    end

    if ~isempty(opts.cmap)
        colormap(ax(i),opts.cmap)
    end
    if ~do3D
    switch opts.draworder
        case 'value'
            hs(i).ZData=cvals(thisgrp); %high expr on top
        case 'random'
            hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
        case 'flat' %not ordered - keep as 2D (eg. for alpha)
        otherwise
    end
    end
end

linkaxes(ax)

% common title/colorbar
ht=[];
if ~isempty(opts.title)
    ht=sgtitle(opts.title); %specify the figure
end
if opts.commonCbar
    rectpos=[opts.margins(1:2),1-opts.margins(3:4)-opts.margins(1:2)];
    hcb=makeCB(ax(end),rectpos,opts.cbgap,opts.cbdims,opts.cbLoc);
    hcb.Limits=[min(cvals),max(cvals)];
    if opts.isdiverging
%         cbtick=[min(cvals),0,max(cvals)];
        cbtick=unique(sort([min(hcb.Ticks),max(hcb.Ticks),0]));
    else
%         cbtick=[min(cvals),max(cvals)];
        cbtick=[min(hcb.Ticks),max(hcb.Ticks)];
    end
    hcb.Ticks=sort(cbtick);
    hcb.TickLength=opts.cbdims(2);
end

hsc.fig=fh;
hsc.ax=ax;
hsc.hs=hs;
hsc.cb=hcb;
hsc.ht=ht;

end
function hcb=makeCB(ax,rectpos,gap,dims,placement)
    switch placement
        case 'east'
            cbh=dims(1);
            cbw=dims(2);
            cbx=rectpos(1)+rectpos(3)+gap;
            cby=rectpos(2)+(rectpos(3)-cbh)/2;
        case 'west'
            cbh=dims(1);
            cbw=dims(2);
            cbx=rectpos(1)-gap-cbw;
            cby=rectpos(2)+(rectpos(4)-cbh)/2;
        case 'north'
            cbw=dims(1);
            cbh=dims(2);
            cbx=rectpos(1)+(rectpos(3)-cbw)/2;
            cby=rectpos(2)+rectpos(4)+gap;
        case 'south'
            cbw=dims(1);
            cbh=dims(2);
            cbx=rectpos(1)+(rectpos(3)-cbw)/2;
            cby=rectpos(2)-gap-cbh;
    end
    hcb=colorbar(ax,placement);
    hcb.Position=[cbx,cby,cbw,cbh];
end