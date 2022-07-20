classdef scScatterViewer < handle
    % use scatter_grp and scatter_val to plot features of a single cell
    % dataset. 

    %TODO: this would be easier if scatter_val and scatter_grp were "chart"
    %objects... setters & getters etc, and all their options as options..

    %TODO: make SCDataset - should revolve around that object.

    %hold on to the data for access via UI
    properties
        %data
        scd %needs a substruct that has .coords
        cells
        genes
        expr

        coords %current coordinates being used (2D/3D)
        current_CData %data for point colors
        current_SData
        is_cat %categorical data? use scatter_grp. user can toggle?
        

        %UI
        fh matlab.ui.Figure
        hsc %scatter_grp/val output struct
    end

    methods
        function obj=scScatterViewer(scd, genes, cells, expr, opts)
        arguments
            scd
            genes
            cells
            expr
            opts
        end
        
        if size(cells,1)>nnz(scd.subset)
            cells=cells(scd.subset,:);
            expr=expr(scd.subset,:);
        end

        obj.scd=scd;
        obj.genes=genes;
        obj.cells=cells;
        obj.expr=expr;

        end

        function plot_groups(obj, newgroups)
        end

        function plot_vals(obj, newgroups)
        end

        function keypress(src, event)
            disp(event.KeyData)
        end
    end
end