function hsc=scatter_grp(coords,groupby,splitby,opts,scopts)
arguments
    coords
    groupby=[]
    splitby=[]
    opts.title=[] %label to use as title
    opts.fig=[] %specify figure to plot in (if empty, new figure)
    opts.ax=[] %specify axis to plot in (ignored if multipanel)
    opts.margins=[0.1,0.1,0.1,0.1] %[left, bottom, right, top]
    opts.tilegaps=[0.05,0.05] %multi-subplot gaps [horizontal, vertical]
    opts.panels=[] %specific layout of [nrow, ncol] subplots

    opts.gcols=[]
    opts.draworder {mustBeMember(opts.draworder,["random","index","flat"])} ='flat'
    scopts.?matlab.graphics.chart.primitive.Scatter
end

[nObs,nDim]=size(coords);

% default markersize
% defsz=12-log(nObs);
defsz=3;

%support 3D scatterplots
do3D=false;
soptargs=namedargs2cell(scopts);
if nDim==2
    scatterfun=@(X,C)scatter(X(:,1),X(:,2),defsz,C,'filled',soptargs{:});
elseif nDim==3
    scatterfun=@(X,C)scatter3(X(:,1),X(:,2),X(:,3),defsz,C,'filled',soptargs{:});
    do3D=true;
%     disp('3D scatterplot')
else
    error('coors are not 2 or 3 dimensional');
end

% organize grouping variable for colors
if isempty(groupby)
    groupby=ones(nObs,1);
end
if ~iscategorical(groupby)
    groupby=categorical(groupby);
end
gcounts=countcats(groupby);
groupby=removecats(groupby);
gnames=categories(groupby);
nGrp=length(gnames);

gcols=opts.gcols;
if size(gcols,1)==length(gcounts)
    gcols=gcols(gcounts>0,:);
else
    if nGrp==1
        gcols=[0.6,0.6,0.6];
    else
        gcols=turbo(nGrp);
    end
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
    ax=opts.ax;
    fh=gcf;
elseif ~isempty(opts.fig)&&isempty(opts.ax) %make new axes
    fh=figure(opts.fig);
    doNewAx=true;
else %both fig and ax specified: splitby only possible if correct # ax
    if length(opts.ax)~=nSplit
        error("Cannot fit "+string(nSplit)+" scatterplots into "+string(length(opts.ax))+" axes");
    end
    fh=figure(opts.fig);
    ax=opts.ax;
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
if length(ax)>nSplit
    delete(ax(nSplit+1:end))
    ax=ax(1:nSplit);
end

for i=1:nSplit
    axes(ax(i))
    hold on
    for j=1:nGrp
        thisgrp=groupby(:)==gnames{j} & splitby(:)==snames{i};
        hs(j)=scatterfun(coords(thisgrp,:),gcols(j,:));
        hs(j).DisplayName=gnames{j};

        if ~do3D
            switch opts.draworder
                case 'index'
                    hs(j).ZData=j*ones(size(hs(j).XData));  %order by category index
                case 'random'
                    hs(j).ZData=rand(size(hs(j).XData));  %randomize the "depth" of points
                case 'flat' %not ordered - keep as 2D (eg. for alpha)
                otherwise
            end
        end
    end
    hold off

    axis off
    if nSplit>1
        title(snames{i})
    end

end

% common title/colorbar
ht=[];
if ~isempty(opts.title)
    ht=sgtitle(opts.title); %specify the figure
end

% ax(end).Selected='on';

if do3D
    axis(ax,'vis3d')
elseif length(ax)>1
    linkaxes(ax)
end


hsc.fig=fh;
hsc.ax=ax;
hsc.hs=hs;
hsc.ht=ht;
hsc.opts=opts;
hsc.scopts=scopts;
end