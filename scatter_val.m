function hsc=scatter_val(coords,cvals,splitby,opts,scopts)
arguments
    coords
    cvals=[]
    splitby=[]
    opts.splitnames=[]
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
    opts.cbJust {mustBeMember(opts.cbJust,["low","mid","high"])}='mid'
    opts.commonCbar=true %true - one cb for all axes same color scale. false - individual cbs
    opts.cbDigits=1

    opts.cmap=[]
    opts.isdiverging=false

    opts.draworder {mustBeMember(opts.draworder,["random","value","valrev","flat"])} ='flat'
    scopts.?matlab.graphics.chart.primitive.Scatter
end
%returns a struct containing all handles for post-modification

[nObs,nDim]=size(coords);

%default colors
if isempty(cvals)
    cvals=ones(nObs,1);
end
[m,n]=size(cvals);
if m~=nObs 
    if n~=nObs
        error("cvals and coords have different number of observations")
    else
        cvals=cvals';
        [m,n]=size(cvals);
    end
end
if n>1
    if ~isempty(splitby)
        warning("cvals represents multiple features: ignoring splitby")
    end
    splitby=categorical(string(1:n));
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
if ~isempty(opts.splitnames) && length(opts.splitnames)==length(snames)
    snames=opts.splitnames;
    splitby=renamecats(splitby,snames);
end
nSplit=length(snames);

% default markersize
% defsz=13-log(nObs);
defsz=3;

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


%decide whether to make figure and/or axes
doNewAx=false;
if isempty(opts.fig)&&isempty(opts.ax) %make new fig+axes
    fh=figure();
    doNewAx=true;
elseif isempty(opts.fig)&&~isempty(opts.ax) %create correct # ax
    if length(opts.ax)~=nSplit
        error("Cannot fit "+string(nSplit)+" scatterplots into "+string(length(opts.ax))+" axes");
    end
    ax=opts.ax(1);
    fh=gcf;
elseif ~isempty(opts.fig)&&isempty(opts.ax) %make new axes
    fh=figure(opts.fig);
    doNewAx=true;
else %both fig and ax specified: splitby only possible if correct # ax
    if length(opts.ax)~=nSplit
        warning("Cannot fit "+string(nSplit)+" scatterplots into "+string(length(opts.ax))+" axes");
    end
    doNewAx=true;
    fh=figure(opts.fig);
end

if doNewAx
    if ~isempty(opts.panels)
        nr=opts.panels(1);
        nc=opts.panels(2);
    else
        nr=floor(sqrt(nSplit));
        nc=ceil(sqrt(nSplit));
    end
    ax=tight_subplot(nr,nc,[],opts.tilegaps,opts.margins([2,4]),opts.margins([1,3]));
end


for i=1:nSplit
    axes(ax(i))

    if n>1
        thisgrp=splitby==snames{i};
        thiscvals=cvals(:,thisgrp);
        thiscoords=coords;
    else
        thisgrp=splitby==snames{i};
        thiscvals=cvals(thisgrp);
        thiscoords=coords(thisgrp,:);
    end
    hs(i)=scatterfun(thiscoords,thiscvals);
    
    if ~opts.commonCbar
        rectpos=ax(i).Position;
        hcb(i)=makeCB(ax(i),rectpos,opts.cbgap,opts.cbdims,opts.cbLoc,opts.cbJust);
        if opts.isdiverging
            ax(i).CLim=abs(max(thiscvals))*[-1,1];
            cbtick=unique([min(thiscvals),0,max(thiscvals)]);
%             cbtick=unique(sort([min(hcb(i).Ticks),max(hcb(i).Ticks),0]));
        else
            ax(i).CLim=[min(thiscvals),max(thiscvals)];
            cbtick=[min(thiscvals),max(thiscvals)];
%             cbtick=[min(hcb(i).Ticks),max(hcb(i).Ticks)];
        end
        hcb(i).Limits=[min(thiscvals),max(thiscvals)];
        cbtick=[ceil(min(cbtick)*10^opts.cbDigits),floor(max(cbtick)*10^opts.cbDigits)]/10^opts.cbDigits;
        hcb(i).Ticks=sort(cbtick);
%         hcb(i).TickLength=opts.cbdims(2);
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
            hs(i).ZData=thiscvals; %high on top
        case 'valrev'
            hs(i).ZData=-thiscvals; %low on top
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
%     if nSplit>1
        ht=sgtitle(opts.title); %specify the figure?
%     else
%         ht=title(opts.title); %specify the axis?
%     end
end
if opts.commonCbar
    rectpos=[opts.margins(1:2),1-opts.margins(3:4)-opts.margins(1:2)];
    hcb=makeCB(ax(end),rectpos,opts.cbgap,opts.cbdims,opts.cbLoc,opts.cbJust);
    hcb.Limits=[min(cvals),max(cvals)];
    if opts.isdiverging
        cbtick=[min(cvals),0,max(cvals)];
%         cbtick=unique(sort([min(hcb.Ticks),max(hcb.Ticks),0]));
    else
        cbtick=[min(cvals),max(cvals)];
%         cbtick=[min(hcb.Ticks),max(hcb.Ticks)];
    end
    cbtick=[ceil(min(cbtick)*10^opts.cbDigits),floor(max(cbtick)*10^opts.cbDigits)]/10^opts.cbDigits;
    hcb.Ticks=sort(cbtick);
    hcb.TickLength=opts.cbdims(2);
end

hsc.fig=fh;
hsc.ax=ax;
hsc.hs=hs;
hsc.cb=hcb;
hsc.ht=ht;
hsc.opts=opts;
hsc.scopts=scopts;

end

%more location options (side,low/center/high)
function hcb=makeCB(ax,rectpos,gap,dims,placement,justify)
    %common features of east/west or north/south
    switch placement
        case {'east','west'}
            cbh=dims(1)*rectpos(4);
            cbw=dims(2);
            switch justify
                case 'low'
                    cby=rectpos(2);
                case 'mid'
                    cby=rectpos(2)+(rectpos(4)-cbh)/2;
                case 'high'
                    cby=rectpos(2)+rectpos(4)-cbh;
            end
        case {'north','south'}
            cbw=dims(1)*rectpos(3);
            cbh=dims(2);
            switch justify
                case 'low'
                    cbx=rectpos(1);
                case 'mid'
                    cbx=rectpos(1)+(rectpos(3)-cbw)/2;
                case 'high'
                    cbx=rectpos(1)+rectpos(3)-cbw;
            end
    end
    %unique features of each side: gap from axis
    switch placement
        case 'east'
            cbx=rectpos(1)+rectpos(3)+gap;
        case 'west'
            cbx=rectpos(1)-gap-cbw;
        case 'north'
            cby=rectpos(2)+rectpos(4)+gap;
        case 'south'
            cby=rectpos(2)-gap-cbh;
    end
    

%     [cbx,cby,cbw,cbh]
    hcb=colorbar(ax,placement);
    hcb.Position=[cbx,cby,cbw,cbh];
end