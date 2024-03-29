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
                    ct.subtypes(i)=[]; %actually delete it
                    ct.nSubtypes=numel(ct.subtypes);
                    return
                end
                if ~isempty(ct.subtypes(i).subtypes)
                    ct.subtypes(i).removeSubtype(name);
                end
            end
            disp(['No CellType with name ' name ' found.'])
        end
        
        function [classID,IDs,tf,markerScore]=classifyByScore(ct,tcounts,genes,method)
            % classify cells into leaf-types of the CellType tree
            %
            % classIDs are in {unclassified, ambiguous, ct.Names()}
            % tf - true if cell satisfies cell type condition, rows=leaf nodes
            %
            % intermediate nodes (i.e. those with children) can have conditions, all their children inherit this
            % condition
            
            IDs=cell(1,size(tcounts,2));
            classID=strings(1,size(tcounts,2));
            this_condition=true(1,size(tcounts,2));
            
            %check this node
            if ~(ct.markers=="")
                if method=="thr"
                    [this_condition, markerScore]=ct.evaluateMarkers(tcounts, genes);
                elseif method=="mean"
                    [this_condition, markerScore]=ct.evaluateMarkers2(tcounts, genes);
                end
            end
            
            if isempty(ct.subtypes)
                tf=this_condition;
                IDs(this_condition)={ct.name};
                return
            end
            
            %otherwise:
            markerScore=[]; %what is desired behavior? return only leaf TF for now. traversing the tree means flattening to leaves... add function that individually evaluates a ct node if desired.
            tf=logical.empty(); 
            tf_parent=this_condition;
            
            %classify subtypes. leaf nodes inherit parent node's markers
            for i=1:ct.nSubtypes
                [~,classID_sub,tf_sub,mscore_sub]=ct.subtypes(i).classifyByScore(tcounts,genes,method);
                tf_sub=repmat(tf_parent,size(tf_sub,1),1) & tf_sub; %parent AND subtype conditions pass
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
        
        
        function [this_condition,markerScore] = evaluateMarkers2(ct, tcounts, genes)
            markerScore=score_genes(ct.markers,tcounts,genes.name,ctrl_size=25); %# per cell
            ct.threshold=geneThresholdOtsu(markerScore);
            this_condition=markerScore>ct.threshold;
        end
        
        function [this_condition,markerScore] = evaluateMarkers(ct, tcounts, genes)
            [gix,ct.markers]=getGeneIndices(ct.markers,genes.name); %note this erases unfound markers!
            markersAbove=tcounts(gix,:)>genes.thr(gix); %shows which genes are above
            markerScore=sum(markersAbove,1); %# per cell
%             scorethr=geneThresholdOtsu(markerScore) %update scorethr?
            this_condition=markerScore>=ct.threshold;
        end
        
        
        function [classID,IDs,tf,markerScore]=iterativeClassifyByScore(ct,tcounts,genes,params,coords)
            %initial partition
            [cellID,IDs,TF,mscores]=ct.classifyByScore(tcounts,genes);

            if params.doImpute
                new_cellID=cellID;
                cellsToImpute=new_cellID=="Unc";
                [new_cellID_imp, newID]=imputeCellType(new_cellID,cellsToImpute,params.impute,coords);
                ix=find(cellsToImpute); 
                new_cellID_imp(ix(newID=="Amb"))="Unc"; %don't impute to "Ambiguous" class.
                numImputed=nnz(new_cellID=="Unc")-nnz(new_cellID_imp=="Unc")
                new_cellID=new_cellID_imp;
                cellID=new_cellID;
            end
            summary(cellID)

            iter=1;
            new_ct=ct;
            while iter<=params.maxIter

                % for each type in cellID that is not Amb/Unc, try to refine markers:
                cats=string(categories(cellID));
                cats(cats=="Amb")=[];
                cats(cats=="Unc")=[];
                for i = 1:length(cats)
                    initial_markers = new_ct.subtypes(i).markers;
                    %check self % and % in every other type for all genes.
                    self=cellID==cats(i);
                    prct_self=sum(tcounts(:,self)>genes.thr,2)./nnz(self)*100;
                    keep=prct_self>params.selfPrctThr;

                    others=setxor(cats,cats(i),'stable');
                    for j=1:length(others)
                        other=cellID==others(j);
                        prct_other=sum(tcounts(:,other)>genes.thr,2)./nnz(other)*100;
                        keep=keep & prct_other<params.otherPrctThr;
                    end

                    markers=genes.name(keep);
                    if params.keep_initial_ids
                        markers=unique([markers(:);initial_markers(:)]);
                    end

                    %set the new markers
                    [gix,markers]=getGeneIndices(markers,genes.name);
                    new_ct.subtypes(i).markers=markers;

                    %set a new score threshold
                    geneThresholds=genes.thr(gix);
                    genesAbove=tcounts(gix,:)>=geneThresholds;
                    markerScore=sum(genesAbove,1);
                    scoreThreshold=round(geneThresholdOtsu(markerScore));
                    scoreThreshold=max(ct.subtypes(i).threshold,scoreThreshold); %don't go below init thr
                    new_ct.subtypes(i).threshold=scoreThreshold;
                end

                % classify based on new markers
                [new_cellID,IDs,TF,mscores]=new_ct.classifyByScore(tcounts,genes);

                if params.doImpute
                    cellsToImpute=new_cellID=="Unc";
                    [new_cellID_imp, newID]=imputeCellType(new_cellID,cellsToImpute,params.impute,coords);
                    ix=find(cellsToImpute); 
                    new_cellID_imp(ix(newID=="Amb"))="Unc"; %don't impute to "Ambiguous" class.
                    numImputed=nnz(new_cellID=="Unc")-nnz(new_cellID_imp=="Unc")
                    new_cellID=new_cellID_imp;
                end
                summary(new_cellID)

                if isequal(cellID,new_cellID)
                    cellID=new_cellID;
                    break
                else
                    cellID=new_cellID;
                end

                iter=iter+1;
            end
        end
        
        
        function [classID,IDs,tf,markerScore,cellClassID]=classifyClusterByScore(ct,clusterID,tcounts,genes)
            % assign all cells of a cluster to the cell type most common in
            % the cluster. (mode)
            [cellClassID,IDs,tf,markerScore]=ct.classifyByScore(tcounts,genes);
            classID = cellClassID;
            ids = unique(clusterID);
            for id = ids
                thisclass = mode(cellClassID(clusterID==id));
                classID(clusterID==id) = thisclass;
            end
        end
        
        function [classID,dom_markers]=classifyClusterByDOM(ct,clusterID,ncounts,tcounts,genes)
            ctnames=ct.Names();
            clusterID=categorical(clusterID);
            clust_names=categories(clusterID);
            
            if nargin-1==2
                %DOM analysis passed in as ncounts
                DOM=ncounts;
                for i=1:length(clust_names)
                    domgenes=DOM.(clust_names{i}).name;
                    for j=1:length(ctnames)
                        dom_markers{i,j}=domgenes(ismember(domgenes,ct.Markers(ctnames{j})));
                    end
                end
            else
                % check dominant genes (no stats) for enrichment of ct.Markers 
                params.direction='up';
                params.prctMetric='diff';
                params.prctESThr=20;
                params.minPrctEffect=0;
                params.fcExprThr=2;
                params.minFcExpr=1.0;
                params.combine_prct_expr='or';

                self.names={'self'};
                self.minprct=20;
                % self.maxprct=15; %down
                self.poolmethod='all';

                other.names={'other'};
                % other.minprct=20; %down
                other.maxprct=5;
                other.poolmethod='some';
                other.some_N=2;

                classID=strings(1,size(tcounts,2));
                genesByClust=getExpression(genes,ncounts,tcounts,clusterID);
                for i=1:length(clust_names)
                    s=self; o=other;
                    s.names=clust_names(i);
                    o.names=setdiff(clust_names, s.names, 'stable');
                    [dominant, ~, ~]=identifyDominantGenes(genesByClust, s, o, [], params);
                    Dnames{i}=dominant;
                    for j=1:length(ctnames)
                        dom_markers{i,j}=dominant.name(ismember(dominant.name,ct.Markers(ctnames{j})));
                    end
                end
            end
            
            numdom=cellfun(@length,dom_markers);
            [nmax,maxix]=max(numdom,[],2);
            for i=1:length(clust_names)
                classID(clusterID==clust_names{i})=ctnames(maxix(i));
            end
            classID=categorical(classID,ctnames,ctnames);
            
            dom_markers=array2table(dom_markers,'RowNames',clust_names,'VariableNames',ctnames);
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
        
        function ct=setMarkers(ct,name,newMarkers,append)
            if append
                markerList=ct.preorderQuery('markers',name);
                newMarkers=unique([markerList(:);newMarkers(:)]);
            end
            ct.preorderSet(name,'markers',newMarkers);
        end
        
        function markerList=Markers(ct,option)
            if ~exist('option','var')
%                 option="leaves";
                option=[];
            end
            markerList=ct.preorderQuery('markers',option);
            markerList(markerList=="")=[];
        end
        
        function ct=setColors(ct,names,colors)
            %names could be 'all' or 'leaves' or list of specific names?
        end
        
        function colorArray=Thresholds(ct,option)
            if ~exist('option','var')
%                 option="leaves";
                option=[];
            end
            colorArray=ct.preorderQuery('threshold',option);
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
            elseif any(strcmp(option,ct.name)) %strcmp to handle empty option
                result=ct.(query);
            elseif any(strcmp(option,"all")) && ~any(ct.markers=="")
                result=ct.(query);
            end
            
            for i=1:ct.nSubtypes
                result=[result;ct.subtypes(i).preorderQuery(query,option)];
            end
        end
        
        function ct=preorderSet(ct, name, prop, val)
            if ct.name==name
                ct.(prop)=val;
            else
                for i=1:ct.nSubtypes
                    ct.subtypes(i)=ct.subtypes(i).preorderSet(name, prop, val);
                end
            end
        end
    end
    
end