function [v,ax,ylh,labs]=violinmatrix(data,group,groupColors,varNames,gap,margins)
% data has rows=observations, columns=variables

if ~iscategorical(group)
    group=categorical(group);
end
groupnames=categories(group);
nGroups=length(groupnames);
nVars=length(varNames);

if isempty(groupColors)
    groupColors=cbrewer('qual','Set1',nGroups);
end

ylh=[];
[~,pos]=tight_subplot(nVars,1,gap,margins(1:2),margins(3:4));
% POS=cat(1,pos{:});
for i=1:nVars
    ax(i)=subplot('Position',pos{i});
    ax(i).YTickLabelMode='auto';
    
%     ax(i)=subplot_tight(nVars,1,i,spmargins);
    v(i,:)=violinplot(data(:,i),group,'BandWidth',0.1);
%     ylh(i)=ylabel(varNames{i});
    
%     if mod(i,2)==0
%         ax(i).YAxisLocation='right';
%     end

    xlim(ax(i),[0.5,nGroups+0.5]);
    ax(i).YAxisLocation='right';
    
    for j=1:nGroups
        v(i,j).ViolinColor=groupColors(j,:);
        v(i,j).ViolinAlpha=0.8;
        v(i,j).ShowData=false;
        v(i,j).MedianPlot.Visible='off';
        v(i,j).BoxPlot.Visible='off';
        v(i,j).WhiskerPlot.Visible='off';
        v(i,j).MedianPlot.Visible='off';
    end
    
%     grid on
    if i<nVars
        ax(i).XAxis.Visible='off';
    end
    
    XLIM=xlim();
    YLIM=ylim();
    labs(i)=text(XLIM(1)-0.05,YLIM(1)+diff(YLIM)/2,varNames(i),'HorizontalAlignment','right','VerticalAlignment','middle');
end
% labs=text(POS(:,1),POS(:,2),varNames,'Units','Normalized','HorizontalAlignment','right','VerticalAlignment','top');
% linkaxes(ax,'y')