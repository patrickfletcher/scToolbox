function [ax, hs, hc]=plotScatter(X,colorby,group,colors,figID,subplotdims,sp_params,docolorbar)
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
marker='o';

nObs=size(X,1);
markerSize=11-log(nObs);

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

if ~exist('docolorbar','var')||isempty(docolorbar)
    docolorbar=true;
end


switch lower(colorby)
    
    case 'group'
        
        if isempty(group)
            group=ones(size(X,1),1);
        end
        
        if ~iscategorical(group)
            group=categorical(group);
        end
        groupNames=categories(group);
        groupCounts=countcats(group);
        
        if ~exist('colors','var')||isempty(colors)||size(colors,1)<length(groupNames)
            colors=cbrewer('qual','Set1',max(length(groupNames),3));
        end
        
        colors=colors(groupCounts>0,:);
        group=removecats(group);
        groupNames=categories(group);
            
        if ~exist('subplotdims','var')||isempty(subplotdims)||~(class(subplotdims)=="matlab.graphics.axis.Axes")
%             ax=subplot(1,1,1);
            ax=tight_subplot(1,1,1,sp_params.gap,sp_params.marg_h,sp_params.marg_w);
        else
            ax=subplotdims;
            axes(ax); %bring specified axes into focus
        end
        
%         hs=gscatter(X(:,1),X(:,2),group);
%         for i=1:length(hs)
%             hs(i).Marker='o';
%             hs(i).MarkerFaceColor=colors(i,:);
%             hs(i).MarkerEdgeColor=colors(i,:)*0.66; %darker edge of same color
%             hs(i).MarkerSize=4; 
%             hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
%         end

    %alt using basic scatter only
        hold on
        for i=1:length(groupNames)
            thisGroup=group==groupNames{i};
            hs(i)=scatter(X(thisGroup,1),X(thisGroup,2),markerSize,colors(i,:),marker,'filled');
            hs(i).MarkerEdgeColor=colors(i,:)*0.66; %darker edge of same color
            hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
            hs(i).DisplayName=groupNames{i};
        end
        hold off

        ax.SortMethod='depth';
        
        legend off %create user legend outside
%         axis tight
        axis off
        
        hc=[];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'value'
        
%         cmap=flipud(cbrewer('div','RdBu',11)); %test binning?
        cmap=cbrewer('seq','OrRd',17); %test binning?
        cmap=[[.7,.85,0.9];cmap(4:end,:)]; %add a teal color for zero; skip lower couple colors

        %group - simply labels for the plots...
        %colors - size(X,1) by number of variables to plot
        
        [m,n]=size(colors);
        if m~=nObs && n~=nObs
            error('Number of observations does not match number of scatter points')
        elseif n==size(X,1) 
            colors=colors'; %force columns
            [~,n]=size(colors);
        end
        nVars=n;
        
        if isempty(group)
            group=repmat("",n,1);
        end
        
        mustMakeAxes=true;
        if ~exist('subplotdims','var')||isempty(subplotdims)
            p=numSubplots(nVars);
            nr=p(1);
            nc=p(2);
        elseif class(subplotdims)=="matlab.graphics.axis.Axes"
            mustMakeAxes=false;
            ax=subplotdims; %could be an array of axes
        elseif isnumeric(subplotdims)
            if prod(subplotdims)>=nVars
                nr=subplotdims(1);
                nc=subplotdims(2);
            else
                p=numSubplots(nVars);
                nr=p(1);
                nc=p(2);
            end
        end
        if mustMakeAxes
            ax=tight_subplot(nr,nc,[],sp_params.gap,sp_params.marg_h,sp_params.marg_w);
        end
        
        for i=1:nVars
            axes(ax(i));
            
            hs(i)=scatter(X(:,1),X(:,2),markerSize,colors(:,i),marker,'filled');
            hs(i).MarkerEdgeColor=hs(i).MarkerFaceColor;
            
            axis equal
            axis tight
            axis off
            
            colormap(cmap);
            
            %custom colobar, small & centered to the right. shows only max/min color
            if docolorbar
                c=colorbar;
                c.Position(1)=ax(i).Position(1)+ax(i).Position(3)+cb_gap;
                c.Position(3)=cb_width;
                c.Position(2)=c.Position(2)+c.Position(4)/4;
                c.Position(4)=c.Position(4)/2;

                c.Ticks=[ceil(10*c.Limits(1)),floor(10*c.Limits(2))]/10; %round to 1 decimal point

                hc(i)=c;
            else
                hc=[];
            end
            
            hs(i).ZData=rand(size(hs(i).XData));  %randomize the "depth" of points
%             hs(i).ZData=colors(:,i); %high expr on top
            
            ax(i).SortMethod='depth';
            
            title(group(i));
            
        end
        
end

if nargout==0
    clear ax
end

end


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
