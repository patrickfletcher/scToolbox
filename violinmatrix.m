function [v,ax,ylh]=violinmatrix(data,group,groupColors,varNames,bandwidth,widths)

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

if ~exist('bandwidth','var')
    bandwidth=0.1;
end
if length(bandwidth)==1
    bandwidth=repmat(bandwidth,1,nVars);
end
% 
% if ~exist('sppars','var')
%     sppars.gap=0.05;
%     sppars.marg_h=0.05;
%     sppars.marg_w=0.05;
% end

if ~exist('widths','var')
    widths=0.3*ones(nGroups,1);
end

t=tiledlayout(nVars,1);
t.Padding='compact';
t.TileSpacing='none';
for i=1:nVars
    ax(i)=nexttile(t);
%     ax(i)=tight_subplot(nVars,1,i,sppars.gap,sppars.marg_h,sppars.marg_w);
    
%     if i<nVars
%         ax(i).XAxis.Visible='off';
%     end
    
    v(i,:)=violinplot(data(:,i),group(:),'BandWidth',bandwidth(i),'Widths',widths);
    for j=1:nGroups
        v(i,j).ViolinColor=groupColors(j,:);
        v(i,j).ViolinAlpha=1;
        v(i,j).ShowData=false;
        v(i,j).MedianPlot.Visible='off';
        v(i,j).BoxPlot.Visible='off';
        v(i,j).WhiskerPlot.Visible='off';
    end
    
%     if mod(i,2)==0
%         ax(i).YAxisLocation='right';
%     end
    ylh(i)=ylabel(varNames{i});
    axis tight
    ax(i).XTickLabel=[];
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
ax(end).YTickLabelMode='auto';