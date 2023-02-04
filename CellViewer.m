classdef CellViewer < matlab.graphics.chartcontainer.ChartContainer & ...
                           matlab.graphics.chartcontainer.mixin.Legend & ...
                           matlab.graphics.chartcontainer.mixin.Colorbar

    % Simple scatter wrapper UI for exploring scRNAseq data
    % - keypress interaction, modal UIfigs for selections

    %hold on to the data for access via UI
    properties
        %data
        coords
        cells
        genes
        expr

        %selection for plotting
        xvar
        yvar
        zvar
        cvar

        cIsCat = false
        group
        groupNames
        groupCols

        cmap = turbo
    end

    properties(Access = private,Transient,NonCopyable)
        Scatter matlab.graphics.chart.primitive.Scatter
    end

    properties(Access = private)
        XData = []
        YData = []
        ZData = []
        CData = []
        TitleText = []
    end

    methods
        function obj = CellViewer(coords, cells, genes, expr, opts)
            arguments
                coords %table of possible coordinates????
                cells %cell metadata
                genes %row-labels of expr matrix
                expr %genes x cells expression matrix
                opts.xvar = []
                opts.yvar = []
                opts.zvar = []
                opts.cvar = []
                opts.group = []
                opts.groupNames = {''}
                opts.groupCols = [0.8,0.8,0.8]
            end

            % Convert x, y, and margin into name-value pairs
            args = {'coords', coords, 'cells', cells, 'genes', genes, 'expr', expr};
            
            % Combine args with user-provided name-value pairs.
            args = [args namedargs2cell(opts)];
            
            % Call superclass constructor method
            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
        end
    end

    methods(Access = protected)
        function setup(obj)

            % Get the axes
            ax = getAxes(obj);
            fh = obj.Parent;
            fh.WindowKeyPressFcn=@obj.keypress;
            
            % Create surface and image objects
            obj.Scatter = scatter(ax,[],[]);
            
            % Configure axes, make colorbar visible
            obj.ColorbarVisible = 'off';
            obj.LegendVisible = 'off';
            axis(ax,'padded')
        end

        function update(obj)
            % Update Data and Colormap
            ax = getAxes(obj);

            validVars = obj.resolvePlotVars();
            if ~validVars
                return
            end

            %do the plotting
            for i=1:length(obj.groupNames)
                obj.Scatter.XData = obj.XData;
                obj.Scatter.YData = obj.YData;
                obj.Scatter.ZData = obj.ZData;
                obj.Scatter.CData = obj.CData;
                colormap(ax,obj.cmap)
            end

            if obj.CIsCat
                obj.ColorbarVisible = 'off';
                obj.LegendVisible = 'on';
            else
                obj.ColorbarVisible = 'on';
                obj.LegendVisible = 'off';
            end
            title(getAxes(obj),obj.TitleText);

        end
    end

    methods

        % make sure the x,y,z,c data plus any groupings are ready for
        % update
        function resolvePlotVars(obj)

%             obj.XData = coords(:,1);
%             obj.YData = coords(:,2);
%             if size(coords,2)>2
%                 obj.ZData = coords(:,3);
%             end
%             obj.CData = ones(size(coords,1));
        end

        % keypress to interactively choose plot vars from tables
        function keypress(obj, src, event)
            disp(event.KeyData)
        end

        %usual figure stuff
        function title(obj,txt) 
            obj.TitleText = txt;
        end
    end
end