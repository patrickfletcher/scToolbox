function [AX, HS, HC, HT]=plotScatter(coords,colorby,groups,colors,figID,subplotdims,sp_params,draworder,docolorbar)

%scatter plot of points in X, colored by either category or values
% X - data, rows are points, colums dimension (2 or 3 only)
% colorby - {'group','value'}.
% ident - if 'group', group (categorical) array; if 'value', name(s) of variables (eg. gene names)
% colors - if 'group', one color per group.
%        - if 'value', one color per point (column vector); if more than one column, one subplot per column


%TODO: interface is horrendous. refactor with arguments block

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
    if ~exist('subplotdims','var')||isempty(subplotdims)||~(class(subplotdims)=="matlab.graphics.axis.Axes")
        figID=figure();
    end
else
    figure(figID);
end

if ~exist('draworder','var')||isempty(draworder)
    draworder='value';
%     draworder='random';
end

%handle subplotting...
doRGB=false;
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
    nPlots=length(groups);
    if n==3 && nPlots==1
        doRGB=true;
    end
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
    clf
    AX=tight_subplot(nr,nc,[],sp_params.gap,sp_params.marg_h,sp_params.marg_w);
end

do3D=false;
if size(coords,2)==2
    scatterfun=@(X,C)scatter(X(:,1),X(:,2),markerSize,C,marker,'filled');
elseif size(coords,2)==3
    scatterfun=@(X,C)scatter3(X(:,1),X(:,2),X(:,3),markerSize,C,marker,'filled');
    do3D=true;
else
    error('coors are not 2 or 3 dimensional');
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
                hs(j)=scatterfun(coords(thisGroup,:),color(j,:));
                hs(j).MarkerEdgeColor=color(j,:)*0.8; %darker edge of same color
                
                if ~do3D
                switch draworder
                    case 'index'
                        hs(j).ZData=j*ones(size(hs(j).XData));  %order by category index
                    case 'random'
                        hs(j).ZData=rand(size(hs(j).XData));  %randomize the "depth" of points
                    case 'flat' %not ordered - keep as 2D (eg. for alpha)
                    otherwise
                end
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
%         cmap=cbrewer('seq','OrRd',50); %test binning?
%         cmap=[[.7,.85,0.9];cmap(10:end,:)]; %add a teal color for zero; skip lower couple colors

        cmap=cbrewer('seq','Reds',64);
        % cmap=cbrewer('seq','OrRd',64);
        cmap=[0.8*[1,1,1];cmap];

        if isempty(groups)
            groups=repmat("",n,1);
        end
        
        for i=1:nPlots
            axes(AX(i));
            
            if doRGB
                cols=colors;
            else
                cols=colors(:,i);
            end
            hs(i)=scatterfun(coords,cols);
            hs(i).MarkerEdgeColor=hs(i).MarkerFaceColor;
            
            %custom colobar, small & centered to the right. shows only max/min color
            if ~exist('docolorbar','var')||isempty(docolorbar)
                docolorbar=true;
            end
            if docolorbar && ~doRGB
                c=colorbar;
                c.Position(1)=AX(i).Position(1)+AX(i).Position(3)+cb_gap;
                c.Position(3)=cb_width;
                c.Position(2)=c.Position(2)+c.Position(4)/4;
                c.Position(4)=c.Position(4)/2;

                d=1; 
                c.Ticks=[ceil((10^d)*c.Limits(1)),floor((10^d)*c.Limits(2))]/ (10^d); %round to 1 decimal point

                hc(i)=c;
            else
                hc=[];
            end
            
            if ~do3D
            switch draworder
                case 'value'
                    if doRGB
%                         vals=sum(colors,2);
                        vals=max(colors,[],2);
                    else
                        vals=colors(:,i);
                    end
                    hs(i).ZData=vals; %high expr on top
                case 'random'
                    hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
                case 'flat' %not ordered - keep as 2D (eg. for alpha)
                otherwise
            end
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

if do3D
    axis vis3d
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
