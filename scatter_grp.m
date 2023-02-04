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

    opts.draw_outline=false %when splitby is used draw all points behind
    opts.outline_ident=[]
    opts.outline_alpha=0.9
    opts.outline_col=0.85*[1,1,1]
    opts.outline_size=[]
    opts.gcols=[]
    opts.darken_edges=0
    opts.darken_strength=0.2 %in [0,1]: how far towards black to interp
    opts.draworder {mustBeMember(opts.draworder,["random","index","revind","flat"])} ='flat'
%     opts.do_keypress=false
    opts.textlabs=false

    opts.hide_axis=true

    scopts.?matlab.graphics.chart.primitive.Scatter
end
%returns a struct containing all handles for post-modification

[nObs,nDim]=size(coords);

% default markersize
defsz=12-log(nObs);
if isfield(scopts,"SizeData")
    defsz=scopts.SizeData;
    scopts=rmfield(scopts,'SizeData');
end
if isempty(opts.outline_size)
    opts.outline_size=defsz*3;
end

%support 3D scatterplots
do3D=false;
soptargs=namedargs2cell(scopts);
if nDim==2
    scatterfun=@(X,S,C)scatter(X(:,1),X(:,2),S,C,'filled',soptargs{:});
elseif nDim==3
    scatterfun=@(X,S,C)scatter3(X(:,1),X(:,2),X(:,3),S,C,'filled',soptargs{:});
    do3D=true;
%     disp('3D scatterplot')
else
    error('coors are not 2 or 3 dimensional');
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


% organize grouping variable for colors
if isempty(groupby)
    groupby=ones(nObs,1);
end

% if opts.do_keypress && size(groupby,2)>1
%     allgroupings=groupby;
%     groupby=groupby(:,1);
% end

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
        nc=ceil(nSplit/nr);
%         nc=ceil(sqrt(nSplit));
    end
    ax=tight_subplot(nr,nc,[],opts.tilegaps,opts.margins([2,4]),opts.margins([1,3]));
end
if length(ax)>nSplit
    delete(ax(nSplit+1:end))
    ax=ax(1:nSplit);
end

if ~opts.draw_outline
    hs0=[];
end

ht=[];
for i=1:nSplit
    axes(ax(i))
    if opts.draw_outline
        if isempty(opts.outline_ident)
            out_ix=true(size(coords,1),1);
        else
            out_ix=opts.outline_ident;
        end
        hs0(i)=scatterfun(coords(out_ix,:), opts.outline_size, opts.outline_col);
        hs0(i).Annotation.LegendInformation.IconDisplayStyle='off';
        hs0(i).ZData=(-nGrp-1)*ones(size(hs0(i).XData));
%         hs0.MarkerFaceAlpha=opts.outline_alpha;
%     else
%         hs0=[];
    end
    this_split=splitby(:)==snames{i};
    hold on
    for j=1:nGrp
        thisgrp=groupby(:)==gnames{j};
        cellsub=thisgrp & this_split;
        hs(j)=scatterfun(coords(cellsub,:), defsz, gcols(j,:));
        hs(j).DisplayName=gnames{j};

        if ~do3D
            switch opts.draworder
                case 'index'
                    hs(j).ZData=j*ones(size(hs(j).XData));  %order by category index
                case 'revind'
                    hs(j).ZData=-j*ones(size(hs(j).XData));  %order by reversed category index
                case 'random'
                    hs(j).ZData=rand(size(hs(j).XData));  %randomize the "depth" of points
                case 'flat' %not ordered - keep as 2D (eg. for alpha)
                otherwise
            end
        end
        
        if opts.darken_edges
            edgecol(1,1)=interp1([0,1],[gcols(j,1),0],opts.darken_strength);
            edgecol(1,2)=interp1([0,1],[gcols(j,2),0],opts.darken_strength);
            edgecol(1,3)=interp1([0,1],[gcols(j,3),0],opts.darken_strength);
            hs(j).MarkerFaceColor=gcols(j,:);
            hs(j).MarkerEdgeColor=edgecol;
            hs(j).LineWidth=0.25;
        end

        if opts.textlabs
            tpos=median(coords(cellsub,:));
            if ~do3D
                tpos(3)=1.1;
                switch opts.draworder
                    case 'index'
                        tpos(3)=nGrp+1; 
                    case 'revind'
                    case 'random'
                    case 'flat' 
                    otherwise
                end
            end
            ht(j)=text(tpos(1),tpos(2),tpos(3),gnames{j},HorizontalAlignment="center",VerticalAlignment="middle");
        end
    end

%     if opts.draw_outline
%         hs0(i).SizeData=hs(1).SizeData*3;
%     end

    hold off

    if nSplit>1
        title(snames{i})
    end

    %by default tight/equal?
%     axis(ax(i),'tight','equal')
    ax(i).XTickLabelMode='auto';
    ax(i).YTickLabelMode='auto';
    if opts.hide_axis
        axis(ax(i),'off')
    end
end

% common title/colorbar
htit=[];
if ~isempty(opts.title)
    htit=sgtitle(opts.title); %specify the figure
end

% ax(end).Selected='on';

if do3D
    axis(ax,'vis3d')
end
if length(ax)>1
    linkaxes(ax)
end

axes(ax(1))

hsc.fig=fh;
hsc.ax=ax;
hsc.hs=hs;
hsc.hs0=hs0;
hsc.ht=ht;
hsc.htit=htit;
hsc.opts=opts;
hsc.scopts=scopts;
end