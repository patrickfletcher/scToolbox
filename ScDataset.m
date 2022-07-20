classdef ScDataset < handle & matlab.mixin.Copyable
    % container class for scRNAseq data
    
    %goals:
    % - encapsulate expression matrix, reductions
    % - keep all data in sync: e.g. subsetting
    % - back via data file, only in sense of don't store the counts on save
    
    %ultimately put analysis/plotting functions here too:
    % - simplifies API by not requiring all the data as function args.
        
    %cell/gene data: struct vs map?? Map is like named list. 
    % struct: scd.dimred.my_dr_name.x -> check fieldnames(scd.dimred)
    % map: scd.dimred("my_dr_name").x -> check keys(scd.dimred)
    
    properties
        countsfile
        doSparse
        isPreV3
        
        organism
        mito_prefix
        bm_dataset
        symbol_db
        
        %indices into original counts matrix (for loadCounts)
        orig_cellsub
        orig_genesub

        %map rows=genes, cols=cells: can be a table - array2table(X,'rownames',genes)
        counts
        ncounts
        tcounts
        
        %original items from countsfile
        geneid
        genename
        cellid
        
        %cell data (sync on cell subsetting)
        % metadata table - can add variables to tables at will
        % reductions: map - scd.dimred("umap").coords  (pca, umap, etc)
        % clusterings: map - eg, scd.clust("leiden_r1.0").clusterID
        cells=table
        neighbors
        dimred
        clust
        
        idents %the default cell-grouping variable
        
        %gene-wise data (sync on gene subsetting)
        % metadata and computed items
        % grouped expression summaries: map - scd.ges('type')
        %variable features - hvg.ix must remain valid (make it logical). support multiple? if so: map
        genes=table
        ges
        hvg

        %arbitrary data: map - scd.misc('qcData')
        misc
        
    end
    
    properties (Dependent, SetAccess=private)
        num_cells
        num_genes
    end
    
    
    %%%%%%%%%%%
    % create/load/save
    methods
        
        %Constructor
        function this=ScDataset(countsfile, organism, doSparse, isPreV3)
            arguments
                countsfile
                organism
                doSparse = true
                isPreV3 = false
            end
            
            switch organism
                case 'hsapiens'
                    this.mito_prefix='MT-';
                    this.symbol_db='hgnc';
                case 'mmusculus'
                    this.mito_prefix='Mt-';
                    this.symbol_db='mgi';
                case 'rnorvegicus'
                    this.mito_prefix='Mt-';
                    this.symbol_db='rgd';
            end
            this.organism=organism;
            this.bm_dataset=organism+"_gene_ensembl";
            
            this.doSparse = doSparse;
            this.isPreV3 = isPreV3;
            this.loadCounts(countsfile)
            
            this.orig_cellsub=true(1,size(this.counts,2));
            this.orig_genesub=true(size(this.counts,1),1);
%             this.dimred=containers.Map();
%             this.clust=containers.Map();
%             this.ges=containers.Map();
%             this.hvg=containers.Map();
%             this.misc=containers.Map();
%             this.ct=containers.Map();
            this.neighbors=struct();
            this.dimred=struct();
            this.clust=struct();
            this.ges=struct();
            this.hvg=struct();
            this.misc=struct();
            this.ct=struct();
        end
        
        function this=saveobj(this)
            this=copy(this);
            %don't store the big matrices; reload them using function
            this.counts=[];
            this.ncounts=[];
            this.tcounts=[];
        end
        
        function loadCounts(this, countsfile, doSparse, isPreV3)
            arguments
                this SCDataset
                countsfile = this.countsfile
                doSparse = this.doSparse
                isPreV3 = this.isPreV3
            end
            
            this.countsfile = countsfile;
            [C,gene_name,barcode,gene_id]=load10Xh5matrix(countsfile,doSparse,isPreV3);
            
            if ~isempty(this.orig_cellsub)
                C=C(:,this.orig_cellsub);
                barcode=barcode(this.orig_cellsub);
            end
            N=normalizeCounts(C); %TODO: if new normalization methods are supported, will need params here
            T=log1p(N); 
            %normalization by total counts happens before subsetting genes... or can I do it after?
            if ~isempty(this.orig_genesub)
                C=C(this.orig_genesub,:);
                N=N(this.orig_genesub,:);
                T=T(this.orig_genesub,:);
                gene_name=gene_name(this.orig_genesub);
                gene_id=gene_id(this.orig_genesub);
            end
            this.counts = C;
            this.ncounts = N;
            this.tcounts = T;
            
            this.geneid=gene_id;
            this.genename=gene_name;
            this.cellid=barcode;
            
            if isempty(this.genes) 
                this.initGenes();
            end
            if isempty(this.cells) %else, assumes cells already matches subset
                this.initCells();
            end
            
        end
        
        function initGenes(this)
            gtab = table;
            gtab.id = this.geneid;
            gtab.orig_name = this.genename; 
            gtab.name=make_vars_unique(this.genename);
            gtab.n_cells = sum(this.counts>0,2);
            gtab.n_umi = sum(this.counts,2);

            gtab.Properties.RowNames=gene_id; %allows fast lookups by id...
            this.genes=gtab;
        end
        
        function updateSymbols(this)
            BM = query_biomart(this.bm_dataset,...
                {'gene_biotype',cellstr(this.symbol_db+"_symbol")},'ensembl_gene_id',this.geneid);
            this.genes = [this.genes,BM(:,2:end)];
            this.genename(BM.external_gene_name~="")=BM.external_gene_name(BM.external_gene_name~="");
            this.genename=make_vars_unique(this.genename);
            this.genes.name=this.genename;
        end
        
        function initCells(this)
            ctab=table;
            ctab.id=this.cellid(:);
            ctab.n_genes=sum(this.counts>0,1)';
            ctab.n_umi=sum(this.counts,1)';
            ctab.log1p_n_genes=log10(ctab.n_genes+1);
            ctab.log1p_n_umi=log10(ctab.n_umi+1);

            %fraction mito genes
            mtGenes=contains(gene_name,this.mito_prefix);
            mtCounts=sum(this.counts(mtGenes,:),1);
            fracMito=mtCounts./ctab.n_umi';
            ctab.fracMT=fracMito(:);

            labs=repmat("",nCells,1);
            for i=1:length(sampleinfo.tags)
                has_tag_j = contains(barcode, sampleinfo.tags(i));
                labs(has_tag_j,1)=sampleinfo.labels(i);
            end
            ctab.sample_orig=categorical(labs);

            ctab.sample=ctab.sample_orig;
            ctab.sample(ctab.sample=="AP1")="whole";
            ctab.sample(ctab.sample=="AP2")="whole";
            ctab.sample(ctab.sample=="Anterior1")="anterior";
            ctab.sample(ctab.sample=="Posterior1")="posterior";
            ctab.sample=removecats(ctab.sample);

            samples=categories(ctab.sample);
        end
        
        function addModuleFraction()
        end
        
        function addCellMetadata()
        end
        
        
    end
    
    % utility methods
    % - initial count summary, fraction mitochondial etc (add to cells, genes)
    % - gene name updates (wrap gprof_convert)
    % - norm/transform ? or just call other functions (for now)
    % - cell/gene subsetting and indexing
    % for now simpler wrappers for:
    % - QC stuff (
    % - dimred stuff (findVarGenes, doPCA, sctransform_pca, neighbors, UMAP, clustering)
    % - group expression summary, DEGs
    methods
        
        % subsref
        function [ix, g]=gene_query(this, g)
            
        end
        
        %fetch data -> look in gene names OR cell table var names...

        % filter cells/genes. applies selection to all relevant properties
        function this=subset(this, cellsub, genesub)
            arguments
                this SCDataset
                cellsub = true(1,this.num_cells)
                genesub = true(this.num_genes,1)
            end
            % TODO: should anything beyond logical index be supported?
            
            % if requested, return the result without modifying the original
            if nargout>0
                this = copy(this);
            end
            
            if ~isempty(this.counts)
                this.counts=this.counts(genesub, cellsub);
            end
            if ~isempty(this.ncounts)
                this.ncounts=this.ncounts(genesub, cellsub);
            end
            if ~isempty(this.tcounts)
                this.tcounts=this.tcounts(genesub, cellsub);
            end
            
            if ~isempty(this.cells)
                this.cells=this.cells(cellsub,:);
            end
            if ~isempty(this.genes)
                this.genes=this.genes(genesub,:);
            end
            
            % args are index into current counts. also update global subs
            cellsub_full_ix = find(this.orig_cellsub);
            this.orig_cellsub(cellsub_full_ix(~cellsub))=false;
            
            genesub_full_ix = find(this.orig_genesub);
            this.orig_genesub(genesub_full_ix(~genesub))=false;
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
        
    end
    
end