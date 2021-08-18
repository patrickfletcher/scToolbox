classdef goatools < handle
    %simple wrapper class for some python goatools functionality
    
    properties
        dag %the whole OBO graph datastructure
    end
    
    methods
        %constructor
        function gt = goatools()
            arguments
                
            end
            GODagargs=pyargs('optional_attrs','relationship');
            gt.dag=py.goatools.obo_parser.GODag(GODagargs);
            
        end
    end
    
end