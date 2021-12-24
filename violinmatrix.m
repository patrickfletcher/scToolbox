function [v,ax,ylh]=violinmatrix(data,group,varNames,groupColors,options)
arguments
    data
    group
    varNames
    groupColors
    options.bandwidths=[]
    options.widths=0.3
    options.ymode='left'
    options.ylabmode='vertical'
    options.Padding='compact'
    options.TileSpacing='none'
    options.ax=[]
end

%thr=scalar: one line across all (gene thresh)
%thr=matrix 2x size nGroups: per violin hi/lo?
% thrw=0.33;

if ~iscategorical(group)
    group=categorical(group);
end
groupcounts=countcats(group);
group=removecats(group);
groupnames=categories(group);
nGroups=length(groupnames);
nVars=length(varNames);

if isempty(groupColors)
    groupColors=lines(nGroups);
else
    groupColors(groupcounts==0,:)=[];
end

widths=[];
if length(options.widths)==1
    widths=repmat(options.widths,1,nGroups);
elseif length(options.widths)==nGroups
    widths=options.widths;
end

if isempty(options.ax)
t=tiledlayout(nVars,1);
t.Padding=options.Padding;
t.TileSpacing=options.TileSpacing;
end

for i=1:nVars
    if isempty(options.ax)
        ax(i)=nexttile(t);
    else
        ax(i)=options.ax(i);
    end
    axes(ax(i))
    
%     if i<nVars
%         ax(i).XAxis.Visible='off';
%     end

    bw=[];
    if length(options.bandwidths)==1
        bw=options.bandwidths;
    elseif length(options.bandwidths)>1
        bw=options.bandwidths(i);
    end
    
    if isempty(bw)
        if isempty(widths)
            v(i,:)=violinplot(data(:,i),group(:));
        else
            v(i,:)=violinplot(data(:,i),group(:),'Widths',widths);
        end
    else
        if isempty(widths)
            v(i,:)=violinplot(data(:,i),group(:),'BandWidth',bw);
        else
            v(i,:)=violinplot(data(:,i),group(:),'BandWidth',bw,'Widths',widths);
        end
    end
            
    for j=1:nGroups
        v(i,j).ViolinColor=groupColors(j,:);
        v(i,j).ViolinAlpha=1;
        v(i,j).ShowData=false;
        v(i,j).MedianPlot.Visible='off';
        v(i,j).BoxPlot.Visible='off';
        v(i,j).WhiskerPlot.Visible='off';
    end
    
    ylh(i)=ylabel(ax(i),varNames{i});
    if options.ylabmode=="horizontal"
        ylh(i).Rotation=0;
        ylh(i).HorizontalAlignment='right';
        ylh(i).VerticalAlignment='middle';
    end
    
    if (options.ymode=="lr" && mod(i,2)==0) || options.ymode=="right"
        ax(i).YAxisLocation='right';
        ylh(i).HorizontalAlignment='left';
    end
    axis tight
    ax(i).XTickLabel=[];
    ax(i).YTickLabelMode='auto';
    xlim(ax(i),[0.5,nGroups+0.5]);
    
%     if ~isempty(thr) && isscalar(thr(i,:))
%         line(ax(i),xlim(ax(i)),[1,1]*thr(i),'color',[0.7,0.7,0.7])
%     else
%         for j=1:nGroups
%             line(ax(i),j+thrw*[-1,1],[1,1]*thr(i),'color',[0.7,0.7,0.7])
%         end
%     end
end
ax(end).XTickLabel=groupnames;
% ax(end).YTickLabelMode='auto';