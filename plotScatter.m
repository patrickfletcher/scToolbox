function [AX, HS, HC, HT]=plotScatter(coords,colorby,groups,colors,figID,subplotdims,sp_params,draworder,docolorbar)
%scatter plot of points in X, colored by either category or values
% X - data, rows are points, colums dimension (2 or 3 only)
% colorby - {'group','value'}.
% ident - if 'group', group (categorical) array; if 'value', name(s) of variables (eg. gene names)
% colors - if 'group', one color per group.
%        - if 'value', one color per point (column vector); if more than one column, one subplot per column


%TODO: interface is horrendous. refactor. inputparser?
%TODO: move away from "exist" method of input checking
%TODO: option - add callback function for lasso select points, return their coords/indices to a workspace var
%TODO: mouse-over function: make it display the type name or expr value

%TODO: pass in args to control marker size etc
%TODO: input che`cking... :/
nObs=size(coords,1);
marker='o';
% markerSize=5;
markerSize=13-log(nObs);

%common setup
if ~exist('sp_params','var')||isempty(sp_params)
%     spmargins=[0.01,0.01];
    sp_params.gap=0.1;
    sp_params.marg_h=0.1;
    sp_params.marg_w=0.1;
    sp_params.cb_gap=0.015;
    sp_params.cb_width=0.015;
elseif ~isstruct(sp_params)
    val=sp_params; clear sp_params;
    sp_params.gap=val(1);
    sp_params.marg_h=val(1);
    sp_params.marg_w=val(1);
    sp_params.cb_gap=val(2);
    sp_params.cb_width=val(2);
end
cb_gap=sp_params.cb_gap;
cb_width=sp_params.cb_width;

if ~exist('figID','var')||isempty(figID)
    figID=figure();
else
    figure(figID);
end

if ~exist('draworder','var')||isempty(draworder)
    draworder='value';
%     draworder='random';
end

%handle subplotting...
if colorby=="group"
    %group: categorical/logical/integer matrix, or cell array if >1 plots
    %colors: color matrix, one row per type. (or cell array of these)
    if iscell(groups)
        nPlots=length(groups);
    else
        [m,n]=size(groups); %should be single row/col vector length nObs
        if m~=nObs && n~=nObs
            error('Number of observations does not match number of scatter points')
        end
        nPlots=1;
        groups={groups};
        if ~iscell(colors)
            colors={colors};
        end
    end
elseif colorby=="value"
    %group - simply labels for the plots...
    %colors - size(X,1) by number of variables to plot
    [m,n]=size(colors);
    if m~=nObs && n~=nObs
        error('Number of observations does not match number of scatter points')
    elseif n==size(coords,1) 
        colors=colors'; %force columns
        [~,n]=size(colors);
    end
    nPlots=n;
end

mustMakeAxes=true; 
p=numSubplots(nPlots);
nr=p(1); nc=p(2);    
if exist('subplotdims','var')&&~isempty(subplotdims)
    if class(subplotdims)=="matlab.graphics.axis.Axes"
        AX=subplotdims; %could be an array of axes
        if length(AX)>=nPlots
            mustMakeAxes=false;
        end
    elseif isnumeric(subplotdims)
        if prod(subplotdims)==nPlots
            nr=subplotdims(1);
            nc=subplotdims(2);
        end
    end
end
if mustMakeAxes
    AX=tight_subplot(nr,nc,[],sp_params.gap,sp_params.marg_h,sp_params.marg_w);
end


switch lower(colorby)
    
    case 'group'
        
        for i=1:nPlots
            axes(AX(i));
            
            group=groups{i};
            color=colors{i};
            
            if isempty(group)
                group=ones(size(coords,1),1);
            end

            if ~iscategorical(group)
                group=categorical(group);
            end
            groupNames=categories(group);
            groupCounts=countcats(group);

            autoColors=true;
            if ~isempty(color)
                if size(color,1)==length(groupNames)
                    color=color(groupCounts>0,:);
                    autoColors=false;
                end
            end
            
            group=removecats(group);
            groupNames=categories(group);
            
            if autoColors
                color=cbrewer('qual','Set1',max(length(groupNames),3));
            end
            
            hold on
            for j=1:length(groupNames)
                thisGroup=group==groupNames{j};
                hs(j)=scatter(coords(thisGroup,1),coords(thisGroup,2),markerSize,color(j,:),marker,'filled');
                hs(j).MarkerEdgeColor=color(j,:)*0.66; %darker edge of same color
                
                switch draworder
                    case 'index'
                        hs(j).ZData=j*ones(size(hs(j).XData));  %order by category index
                    otherwise
                        hs(j).ZData=rand(size(hs(j).XData));  %randomize the "depth" of points
                end
                
                hs(j).DisplayName=groupNames{j};
            end
            hold off

%             axis(AX(i),'equal');
%             axis(AX(i),'tight');
            axis(AX(i),'off');
            AX(i).SortMethod='depth';

            HC=[];
            if nPlots>1
                HS{i}=hs;
            else
                HS=hs;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'value'
        
%         cmap=flipud(cbrewer('div','RdBu',11)); %test binning?
        cmap=cbrewer('seq','OrRd',17); %test binning?
        cmap=[[.7,.85,0.9];cmap(4:end,:)]; %add a teal color for zero; skip lower couple colors

        if isempty(groups)
            groups=repmat("",n,1);
        end
        
        for i=1:nPlots
            axes(AX(i));
            
            hs(i)=scatter(coords(:,1),coords(:,2),markerSize,colors(:,i),marker,'filled');
            hs(i).MarkerEdgeColor=hs(i).MarkerFaceColor;
            
            %custom colobar, small & centered to the right. shows only max/min color
            if ~exist('docolorbar','var')||isempty(docolorbar)
                docolorbar=true;
            end
            if docolorbar
                c=colorbar;
                c.Position(1)=AX(i).Position(1)+AX(i).Position(3)+cb_gap;
                c.Position(3)=cb_width;
                c.Position(2)=c.Position(2)+c.Position(4)/4;
                c.Position(4)=c.Position(4)/2;

                d=2; 
                c.Ticks=[ceil((10^d)*c.Limits(1)),floor((10^d)*c.Limits(2))]/ (10^d); %round to 1 decimal point

                hc(i)=c;
            else
                hc=[];
            end
            
            switch draworder
                case 'value'
                    hs(i).ZData=colors(:,i); %high expr on top
                otherwise
                    hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
            end
            
%             axis(AX(i),'equal');
%             axis(AX(i),'tight');
            axis(AX(i),'off');
            AX(i).SortMethod='depth';
            colormap(cmap);
            
            title(groups(i));
            if nPlots>1
                HS=hs;
                HC=hc;
            else
                HS=hs;
                HC=hc;
            end
        end
        
end

if nargout==0
    clear AX
end


% end


% function []=parseArgs(t, X, varargin)
% 
% %default parameters
% defaultF=[0.5,0.4]; %near halfmax
% 
% p=inputParser;
% addRequired(p,'t',@(x) isreal(x));
% addOptional(p,'f',defaultF);  %somehow adding the validation function here messes things up
% addParameter(p,'ThresholdPercentiles',defaultthrPtiles,@(x) isreal(x) && numel(x)==2);
% 
% parse(p,t,X,varargin{:});
% 
% f=p.Results.f;
% 
% end
