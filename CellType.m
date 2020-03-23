classdef CellType %< handle & matlab.mixin.Copyable
    %a class to represent a cell type
    % - subtypes (array of CellTypes) provides support for hierarchy definition
    % - pre-order search the tree for info, or to classify.
    
    properties
        name=""
        markers=""
        threshold=0
        color=[0.5,0.5,0.5]
        subtypes=CellType.empty() %recursively defined tree
        nSubtypes=0
    end
    
    methods
        function ct=CellType(name,markers,threshold,color)
            if exist('name','var') && ~isempty(name)
                ct.name=string(name);
            end
            if exist('markers','var') && ~isempty(markers)
                ct.markers=string(markers(:));
            end
            if exist('threshold','var') && ~isempty(threshold)
                ct.threshold=threshold(:);
            end
            if exist('color','var') && ~isempty(color)
                ct.color=color;
            end
        end
        
        function ct=addSubtype(ct,varargin)
            
            if class(varargin{1})=="CellType"
                newct=varargin{1};
            else
                newct=CellType(varargin{:});
            end
            
            %enforce unique names?
            existingNames=ct.Names();
            if ismember(newct.name,existingNames)
                disp('A cell type with name "'+newct.name+'" already exists. Aborting...');
                return
            end
                
            %TODO: for some reason always creates "ans" in base workspace equal to the new cell type
            ct.subtypes(end+1)=newct;
            ct.nSubtypes=numel(ct.subtypes);
        end
        
        function ct=removeSubtype(ct,name)
            %remove first cell type found with given name, pre-order search
            for i=1:ct.nSubtypes
                if name==ct.subtypes(i).name
                    %delete this subtype, but elevate any children of to-be-deleted node to subtypes of its parent
                    if ~isempty(ct.subtypes(i).subtypes)
                        if i==1
                            ct.subtypes=[ct.subtypes(i).subtypes,ct.subtypes(2:end)]; %add to the beginning
                        elseif i==length(ct.nSubtypes)
                            ct.subtypes=[ct.subtypes(1:end-1),ct.subtypes(i).subtypes]; %add to the end
                        else
                            ct.subtypes=[ct.subtypes(1:i-1),ct.subtypes(i).subtypes,ct.subtypes(i+1:end)];
                        end
                    end
                    
                    ct.nSubtypes=numel(ct.subtypes);
                    return
                end
                if ~isempty(ct.subtypes(i).subtypes)
                    ct.subtypes(i).removeSubtype(name);
                end
            end
            disp(['No CellType with name ' name ' found.'])
        end
        
        function [classID,IDs,tf,markerScore]=classifyByScore(ct,tcounts,genes)
            % classify cells into leaf-types of the CellType tree
            %
            % classIDs are in {unclassified, ambiguous, ct.Names()}
            % tf - true if cell satisfies cell type condition, rows=leaf nodes
            %
            % intermediate nodes (i.e. those with children) can have conditions, all their children inherit this
            % condition
            
            IDs=cell(1,size(tcounts,2));
%             classID=cell(1,size(tcounts,2));
            classID=strings(1,size(tcounts,2));
            
            %keep intermediate nodes separate
            this_condition=true(1,size(tcounts,2));
            if ~(ct.markers=="")
                [gix,ct.markers]=getGeneIndices(ct.markers,genes.name);
                M=tcounts(gix,:);
                markersAbove=M>=genes.thr(gix); %shows which genes are above
                markerScore=sum(markersAbove,1); %# per cell
%                 scoreThreshold=geneThresholdOtsu(markerScore);
                scoreThreshold=ct.threshold;
                this_condition=markerScore>=scoreThreshold;
            end
            
            if isempty(ct.subtypes)
                tf=this_condition;
                IDs(this_condition)={ct.name};
                return
            end
            
            %otherwise:
            tf=logical.empty();
            tf_parent=this_condition;
            markerScore=[];
            
            %classify subtypes. leaf nodes inherit parent node's markers
            for i=1:ct.nSubtypes
                [~,classID_sub,tf_sub,mscore_sub]=ct.subtypes(i).classifyByScore(tcounts,genes);
                tf_sub=repmat(tf_parent,size(tf_sub,1),1)&tf_sub;
                tf=[tf ; tf_sub];
                markerScore=[markerScore; mscore_sub];
                ix=find(any(tf_sub,1));
                for j=ix
                    IDs{j}=[IDs{j},classID_sub{j}];
                end
            end
            
            %resolve ambiguous and unclassified
%             nID=cellfun(@numel,IDs);
            nID=sum(tf,1);
            unc=nID==0;
            uniqueClass=nID==1;
            amb=nID>1;
            
            %return classID as categorical
            classID(uniqueClass)=IDs(uniqueClass);
            if any(amb)
                classID(amb)="Amb";
            end
            if any(unc)
                classID(unc)="Unc";
            end
            catnames=[ct.Names;"Amb";"Unc"]; %to keep hierarchy ordering
            classID=categorical(classID,catnames,catnames);
        end
        
        function [classID,cellClassID,IDs,tf]=classifyClusterByScore(ct,clusterID,tcounts,genes)
            % assign all cells of a cluster to the cell type most common in
            % the cluster. (mode)
            [cellClassID,IDs,tf]=ct.classifyByScore(tcounts,genes);
            classID = cellClassID;
            ids = unique(clusterID);
            for id = ids
                thisclass = mode(cellClassID(clusterID==id));
                classID(clusterID==id) = thisclass;
            end
        end
        
        function ct=setNames(ct,oldNames,newNames)
        end
        
        function nameArray=Names(ct,option)
            if ~exist('option','var')
%                 option="leaves";
                option=[];
            end
            nameArray=ct.preorderQuery('name',option);
        end
        
        function ct=setMarkers(ct,names,newMarkers,newThreshold,append)
        end
        
        function markerList=Markers(ct,option)
            if ~exist('option','var')
%                 option="leaves";
                option=[];
            end
            markerList=ct.preorderQuery('markers',option);
        end
        
        function ct=setColors(ct,names,colors)
            %names could be 'all' or 'leaves' or list of specific names?
        end
        
        function colorArray=Colors(ct,option)
            if ~exist('option','var')
%                 option="leaves";
                option=[];
            end
            colorArray=ct.preorderQuery('color',option);
        end
        

        %TODO: plot - graph plot showing hierarchy - layout layered, root=source, leaves=sinks
        function plot(ct)
        end
        function G=asGraph(ct)
            %determine the graph adjacency matrix - global representation
            
        end
        
    end
    
    methods (Access=private)
        function result=preorderQuery(ct,query,option)
            
            result=[];
            %if has subtypes, return this result only if "all" were requested, otherwise returns only leaves of the tree
            if isempty(option) && isempty(ct.subtypes)
                result=ct.(query);
            elseif any(strcmp(option,ct.name))
                result=ct.(query);
            end
            
            for i=1:ct.nSubtypes
                result=[result;ct.subtypes(i).preorderQuery(query,option)];
            end
        end
    end
    
end