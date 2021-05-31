classdef EnrichPlot < matlab.graphics.chartcontainer.ChartContainer ... 
                      & matlab.graphics.chartcontainer.mixin.Legend ...
                      & matlab.graphics.chartcontainer.mixin.Colorbar
    
    % A chart class for enrichment dot plots.
    % 
    %
    % should inherit from a more general dotplot?
    properties
        XData (1,:) double = NaN
        YData (1,:) double = NaN
        SizeData (1,:) double = NaN
        CData (1,:) double = NaN
        TitleText (:,:) char = ''
        ColorbarTitle (:,:) char = ''
        SizeLegendTitle (:,:) char = ''
    end
    
    properties (Access = private, Transient, NonCopyable)
        BubbleObject (1,1) matlab.graphics.chart.primitive.Bar
        ErrorBarObject (1,1) matlab.graphics.chart.primitive.ErrorBar
    end
    
    properties (Dependent)
        % Provide properties to support setting & getting
        XLimits (1,2) double
        XLimitsMode {mustBeMember(XLimitsMode,{'auto','manual'})}
        YLimits (1,2) double
        YLimitsMode {mustBeMember(YLimitsMode,{'auto','manual'})}
    end
    %inherited:
    % methods
    %     setup
    %     update
    %     getAxes 
    %     getLayout
    %     getLegend
    %     getColorbar
    % end
    

    methods (Access = protected)
        function setup(obj)
            ax = getAxes(obj);
            obj.BarObject = bar(ax,NaN,NaN);
            hold(ax,'on')
            obj.ErrorBarObject = errorbar(ax,NaN,NaN,NaN);
            obj.ErrorBarObject.LineStyle = 'none';
            obj.ErrorBarObject.LineWidth = 2;
            obj.ErrorBarObject.Color = [0.6 0.7 1];
            hold(ax,'off');
        end
        function update(obj)
            % Update Bar and ErrorBar XData and YData
            obj.BarObject.XData = obj.XData;
            obj.BarObject.YData = obj.YData;
            obj.ErrorBarObject.XData = obj.XData;
            obj.ErrorBarObject.YData = obj.YData;
            
            % Update ErrorBar delta values
            obj.ErrorBarObject.YNegativeDelta = obj.EData;
            obj.ErrorBarObject.YPositiveDelta = obj.EData;
            
            % Update axes title
            ax = getAxes(obj);
            title(ax,obj.TitleText);
        end
    end
    
    methods
        %constructor:
        % - EnrichPlot(x,y,siz,col,...)
        function obj = EnrichPlot(x,y,siz,col,varargin)
            arguments
                x (1,:) double = NaN
                y (1,:) double = NaN
                siz (1,:) double = NaN
                col (1,:) double = NaN
            end
            arguments (Repeating)
                varargin
            end

            % Convert x, y, and margin into name-value pairs
            args = {'XData', x, 'YData', y};
            args = [args varargin];

            % Call superclass constructor method
            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
        end
        
        % xlim method
        function varargout = xlim(obj,varargin)
            ax = getAxes(obj);
            [varargout{1:nargout}] = xlim(ax,varargin{:});
        end
        % ylim method
        function varargout = ylim(obj,varargin)
            ax = getAxes(obj);
            [varargout{1:nargout}] = ylim(ax,varargin{:});
        end
        % title method
        function title(obj,txt)
            obj.TitleText = txt;
        end
        
        % set and get methods for XLimits and XLimitsMode
        function set.XLimits(obj,xlm)
            ax = getAxes(obj);
            ax.XLim = xlm;
        end
        function xlm = get.XLimits(obj)
            ax = getAxes(obj);
            xlm = ax.XLim;
        end
        function set.XLimitsMode(obj,xlmmode)
            ax = getAxes(obj);
            ax.XLimMode = xlmmode;
        end
        function xlm = get.XLimitsMode(obj)
            ax = getAxes(obj);
            xlm = ax.XLimMode;
        end
        
        % set and get methods for YLimits and YLimitsMode
        function set.YLimits(obj,ylm)
            ax = getAxes(obj);
            ax.YLim = ylm;
        end
        function ylm = get.YLimits(obj)
            ax = getAxes(obj);
            ylm = ax.YLim;
        end
        function set.YLimitsMode(obj,ylmmode)
            ax = getAxes(obj);
            ax.YLimMode = ylmmode;
        end
        function ylm = get.YLimitsMode(obj)
            ax = getAxes(obj);
            ylm = ax.YLimMode;
        end
        
        function xticks(obj,ticks)
        end
        function yticks(obj,ticks)
        end
        function xticklabels(obj,ticklabels)
        end
        function yticklabels(obj,ticklabels)
        end
    end
    
                      
end