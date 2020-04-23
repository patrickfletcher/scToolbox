classdef upsetplot
    %for set intersection visuals of >3 sets (or alternative to venn for 3)
    
    %make something like: https://github.com/ImSoErgodic/py-upset, https://github.com/hms-dbmi/UpSetR
    
    properties
        %make separate classes containing each plot
        ax_intersect %bar showing size of intersections
        ax_combs %line/dot plot showing set combination for each bar
        ax_setcount %bar showing size of each set
    end
    
    methods 
        function hup=upsetplot()
        end
    end
end