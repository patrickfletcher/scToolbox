classdef upsetplot
    %for set intersection visuals of >3 sets (or alternative to venn for 3)
    
    %make something like: https://github.com/ImSoErgodic/py-upset, https://github.com/hms-dbmi/UpSetR
    
    properties
        ax_intersect %bar showing size of intersections
        h_bar_intersect %the bar chart object
        
        ax_combs %line/dot plot showing set combination for each bar
        h_comb_plot %the set of line plots indicating combinations
        h_comb_bg %gray background line plots with wider marker/lines
        
        ax_setcount %bar showing size of each set
    end
    
    methods 
        function hup=upsetplot()
        end
    end
end