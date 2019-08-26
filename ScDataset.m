classdef ScDataset < handle & matlab.mixin.Copyable
    % container class for scRNAseq data and its analysis
    
    % - recompute size factors for normalization every time genes/cells are subset?
    
    properties
        countsfile
        
        %map rows=genes, cols=cells: can be a table: array2table(X,'rownames',genes)
        counts
        ncounts
        tcounts
        
        numGenes
        numCells
        
        countThreshold=0 %above this is considered detected
        minCellsPerGene=1
        
        %gene table - can add variables to tables at will
        features=table
        
        %cell table
        cells=table
        
        %cell array to store arbitrary metadata
        metadata={};
        
        %dimensionality reductions (includes params used etc)
        hvgix %indices of highly variable genes
        
        %dr
        pca=struct('coords',[]);
        tsne
        
        %clusterings
        clust
        
        %"cellTypeSet" object to contain array of cellTypes and
        %associated data/methods.  
        cellTypeSet=CellTypeSet.empty;
        
        
    end
    
    
    %%%%%%%%%%%
    methods
        
        %Constructor
        function this=ScDataset(countsfile,counts,feature_id,feature_name,cells,countThreshold)
            if exist('countsfile','var')||~isempty(countsfile)
                this.countsfile=countsfile; %support H5, folder with MTX, mat file?
            end
            if exist('counts','var')||~isempty(counts)
                this.counts=counts;
            end
            if exist('gene_id','var')||~isempty(feature_id)
                this.features.id=feature_id;
            end
            if exist('gene_name','var')||~isempty(feature_name)
                this.features.name=feature_name;
            end
            if exist('cells','var')||~isempty(cells)
                this.cells=cells;
            end
            if exist('countThreshold','var')||~isempty(countThreshold)
                this.countThreshold=countThreshold;
            end
        end
        
        
        function this=saveobj(this)
            this=copy(this);
            %don't store the big matrices; reload them from file upon load
            this.counts=[];
            this.ncounts=[];
            this.tcounts=[];
        end
        
        
        function setCounts(this,counts)
        end
        
        function computeCountSums(this)
            countMask=this.counts>0;
            this.cellData.genesPerCell=sum(countMask,1);
            this.cellData.molecPerCell=sum(this.counts,1);
            this.geneData.cellsPerGene=sum(countMask,2);
            this.geneData.molecPerGene=sum(this.counts,2);
        end
        
        
        function normalizeCounts(this)
            %per cell - divide by each cell's total molec count, then
            %multiply by population median molec count
            this.ncounts=this.counts./repmat(this.cellData.molecPerCell,this.numGenes,1).*median(this.cellData.molecPerCell);
        end
        
        
        function transformCounts(this,base,pseudocount)
            %variance stabilizing transform
            
            %log, choose base (default=10) and pseudocount (default=1)
            
            if ~exist('base','var')||isempty(base)
                base=10;
            end
            if ~exist('pseudocount','var')||isempty(pseudocount)
                pseudocount=1;
            end
            
            switch base
                case 2
                    this.tcounts=log2(this.counts+pseudocount);
                case 10
                    this.tcounts=log10(this.counts+pseudocount);
                otherwise
                    error('unsupported base for logarithm')
            end
            
            %Freeman-Tukey
            
        end
        
        
        %%% extract subset of the counts matrix
        function result=subset(this, source, genes, cells)
            %genes can be index, or names/IDs (char/str/cell array of char)
            %cells can be index, or cellID
            %source is raw/norm/lognorm
        end
        
    
    %%%%%%%%%%% Plotting methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cellPlot(this)
        end
        
        function genePlot(this)
        end
        
        function heatmap(this,geneNames)
        end
        
    end
    
    
    
    %%%%%%%%%%%
    methods(Static)
        function lObj=loadobj(obj)
            lObj=obj;
            if ~isempty(obj.filename)
                %try to load data from the file
                %                 lObj.counts=load10Xh5matrix(obj.filename,obj.h5path,obj.doSparse);
            else
                %do nothing?
            end
        end
        
    end
    
end