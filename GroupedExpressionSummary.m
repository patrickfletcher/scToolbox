classdef GroupedExpressionSummary
    %data class to compute and store gene expression data summarized for a
    %grouping variable
       
    % in SCDataset, I can wrap this class with SCD methods 
    %scd would already have: gene info, n/tcounts, cell metadata.  could
    % -> just pass in a grouping variable + options
    
    %DEG component: two use cases. Pairwise tests between groups, and
    %pairwise tests between each group and others pooled.

    % for now don't consider thresholds. Pass tcounts=Tthr if desired
    
    % tables vs matrices? get function to extract custom subsets of
    % groups/stats/etc. into a table via "get" method
    
    properties 
        
        %pooled per-gene values (available without group) - one table?
        id %gene id
        name %gene symbol
        gene_count
        cell_count
        cellsub
        genesub
        pooled table %n_cells, n_umi, stats, is_hvg, ...
        
        %grouping variable
        group %grouping variable (value = missing supported)
        n_groups
        group_names %cols of per-group matrices
        group_counts
        pw_combs    %pairwise combinations of groupnames (cols of pw matrices)
        
        %group-wise quantities - would be nice to just add arbitrary slots
        % [nGenes x nGroups], rows=genes, cols=groups        
        mean_n table   %mean normalized counts
        mean_t table    %mean log normalized counts
%         mean_rank table %this is actual effect size measured by kw/rs
        prct table      %>0
        asinp table      %arcsine transformed fraction expressing
        
        %should also do variance-like measures?
        % std/var, fano, cv, dispersion
        
        %%% statistical tests
        %global tests among groups
        p_anova table %columns - kw: mean ranks, one: anova1, means
        p_chi2 table  %chi2 test for proportions
        
        %pairwise contrast statistics - tables, varname=test name?
        p_mc table %result from multcompare following anova
        p_z table %pairwise z-tests for proportions
        
        
        %pairwise effect sizes - compute all pairwise + lookup, vs use
        %functions: compute for requested self/other groups?
        %functions are pretty fast
%         fc_n
%         fc_t
%         lfc
%         d_prct
%         cohen_h
    end
   
    
    methods
        function ges = GroupedExpressionSummary(genes, counts, ncounts, tcounts, group)
            arguments
                genes table
                counts {mustBeNumeric, mustBeEqualRows(counts,genes)}
                ncounts {mustBeNumeric, mustBeEqualSize(ncounts,counts)}
                tcounts {mustBeNumeric, mustBeEqualSize(tcounts,counts)}
                group {groupLengthCheck(group,counts)}
            end
            %constructor.
            % -genes is a table containing the basics: id, name, n_cells, n_umi
            % -ncounts, tcounts: normalized, log1p transformed counts [ng,nc]
            % -group (optional): grouping variable - length nc vector (can have missing values)
            % -threshold_group + thr: Otsu thresholds must be computed on full raw
            % dataset, they are for ambient RNA...
            
            ges.id=genes.id;
            ges.name=genes.name;
            [ges.gene_count, ges.cell_count]=size(counts);
            
            ges.cellsub=~ismissing(group);
            ges.genesub=true(ges.gene_count,1);
            
            if ~isequal(nnz(ges.cellsub),ges.cell_count)
                counts(:,~ges.cellsub)=[];
                ncounts(:,~ges.cellsub)=[];
                tcounts(:,~ges.cellsub)=[];
                group(~ges.cellsub)=[];
                ges.cell_count=length(group);
            end
            
            %compute group summaries
            ges=ges.groupBy(counts, ncounts, tcounts, group);
        end
        
        function ges=groupBy(ges, counts, ncounts, tcounts, group)
            arguments
                ges GroupedExpressionSummary
                counts {mustBeNumeric, mustBeEqualSizeCounts(counts,ges)} 
                ncounts {mustBeNumeric, mustBeEqualSizeCounts(ncounts,ges)} 
                tcounts {mustBeNumeric, mustBeEqualSizeCounts(tcounts,ges)} 
                group {groupLengthCheck(group,ncounts)}
            end
            % add or modify the grouping variable.
            % populates the group-wise descriptive statistics - no thresh
            
            if ~iscategorical(group)
                group=categorical(group);
            end
            group=removecats(group);
            
            
            ges.group=group;
            ges.group_names=categories(group);
            ges.group_counts=countcats(group);
            ges.n_groups=length(ges.group_names);
            ges.pw_combs=allcomb(ges.group_names);
            
            tic
            %get the pooled stats restricted to valid cells
            ges.pooled.name=ges.name;
            ges.pooled.Properties.RowNames=ges.id;
            ges.pooled.n_cells=sum(counts>0,2);
            ges.pooled.n_umi=sum(counts,2);
            ges.pooled.mean_n=mean(ncounts,2);
            ges.pooled.mean_t=mean(tcounts,2);
            ges.pooled.var_n=var(ncounts,[],2);
            ges.pooled.disp=(ges.pooled.var_n-ges.pooled.mean_n)./(ges.pooled.mean_n).^2;
%             ges.pooled.var_t=var(tcounts,[],2);
            toc
            
            tic
            ges.mean_n.id=ges.id;
            ges.mean_t.id=ges.id;
            ges.prct.id=ges.id;
            ges.asinp.id=ges.id;
            ges.mean_n.name=ges.name;
            ges.mean_t.name=ges.name;
            ges.prct.name=ges.name;
            ges.asinp.name=ges.name;
            for i=1:ges.n_groups
                thisName=ges.group_names{i};
                thisGroup=group==thisName;
                ges.mean_n.(thisName)=mean(ncounts(:,thisGroup),2);
                ges.mean_t.(thisName)=mean(tcounts(:,thisGroup),2);
                frac=sum(tcounts(:,thisGroup)>0,2)./ges.group_counts(i);
                ges.prct.(thisName)=frac*100;
                ges.asinp.(thisName)=asin(sqrt(frac));
            end
            ges.mean_n.Properties.RowNames=ges.id;
            ges.mean_t.Properties.RowNames=ges.id;
            ges.prct.Properties.RowNames=ges.id;
            ges.asinp.Properties.RowNames=ges.id;
            toc
        end
        
        
        function markers = findMarkers(ges, selfnames, othernames, options)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                options.direction = 'up'
                options.selfpool = 'all'
                options.otherpool = 'all'
                options.fc_options = []  %pass a struct if desired
                options.d_prct_options = []
                options.fc_thr = 1
                options.prct_diff = 0
                options.min_self_prct = 20
                options.max_self_prct = 100
                options.min_other_prct = 0
                options.max_other_prct = 5
                options.anova_thr = 1
                options.p_thr = 1
                options.chi2_thr = 1
                options.z_thr = 1
            end
            %apply conditions groupwise, pairwise to get list of genes
            
            %self, other %
            self_prct=ges.prct(:,selfnames);
            other_prct=ges.prct(:,othernames);
            
            
            %get the effect sizes
            if isempty(options.fc_options)
                fc=ges.getFoldChanges(selfnames,othernames);
            else
                fc=ges.getFoldChanges(selfnames,othernames, options.fc_options);
            end
            if isempty(options.d_prct_options)
                d_prct=ges.getPrctDiff(selfnames,othernames);
            else
                d_prct=ges.getPrctDiff(selfnames,othernames, options.d_prct_options);
            end
            
            %effect size selection
            
            
            %get DE p-values
            thisCombs=all(ismember(ges.pw_combs,selfnames)|ismember(ges.pw_combs,othernames),2); %index into columns of the P-matrices
            
        end
        
        
        %return tables? or just matrices?
        function result = getFoldChanges(ges, selfnames, othernames, options)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                options.genes = true(ges.gene_count,1)
                options.slot = 'norm'
                options.transform = ''
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            
            %TODO: need a function that processes multiple types of gene query
            genes = ges.gene_subset(options.genes);
            
            switch whichExpr
                case 'norm' %mean before log
                    switch transform_method
                        case 'log10'
                            %log10( mean(ncounts(self))/mean(ncounts(other)) )
                            fc_fun=@(self,other) log10(ges.mean_n{genes,self}./ges.mean_n{genes,other});
                        case 'log10p1'
                            %log10( (mean(ncounts(self))+1)/(mean(ncounts(other))+1) )
                            fc_fun=@(self,other) log10( (ges.mean_n{genes,self}+1)./(ges.mean_n{genes,other}+1) );
                        case 'log2'
                            %log2( mean(ncounts(self))/mean(ncounts(other)) )
                            fc_fun=@(self,other) log2(ges.mean_n{genes,self}./ges.mean_n{genes,other});
                        case 'log2p1'
                            %log2( (mean(ncounts(self))+1)/(mean(ncounts(other))+1) )
                            fc_fun=@(self,other) log2( (ges.mean_n{genes,self}+1)./(ges.mean_n{genes,other}+1) );
                        case 'log1p'
                            fc_fun=@(self,other) log1p(ges.mean_n{genes,self}./ges.mean_n{genes,other});
                        otherwise
                            %mean(ncounts(self))/mean(ncounts(other))
                            fc_fun=@(self,other) ges.mean_n{genes,self}./ges.mean_n{genes,other};
                    end
                    
                case 'lognorm' %log before mean
                    %mean(log10(ncounts(self)+1))/mean(log10(ncounts(other)+1))
                    fc_fun=@(self,other) ges.mean_t{genes,self}./ges.mean_t{genes,other};
            end
            
            result=table();
            for i=1:length(selfnames)
                this_self=selfnames(i);
                for j=1:length(othernames)
                    this_other=othernames(j);
                    result.(this_self+"_"+this_other)=fc_fun(this_self, this_other);
                end
            end
            
%             result=result{:,:}; %for matrix
        end
        
        
        function result = getCohenH(ges, selfnames, othernames, options)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                options.genes = true(ges.gene_count,1)
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            
            h_fun=@(self,other) 2*(ges.asin{options.genes,self}-ges.asin{options.genes,other});
                    
            result=table();
            for i=1:length(selfnames)
                this_self=selfnames(i);
                for j=1:length(othernames)
                    this_other=othernames(j);
                    result.(this_self+"_"+this_other)=h_fun(this_self,this_other);
                end
            end
            
            result=result{:,:}; %hack for now
        end
        
        
        function result = getPrctDiff(ges, selfnames, othernames, options)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                options.genes = true(ges.gene_count,1)
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            
            d_fun=@(self,other) ges.prct{options.genes,self}-ges.prct{options.genes,other};
            
            result=table();
            for i=1:length(selfnames)
                this_self=selfnames(i);
                for j=1:length(othernames)
                    this_other=othernames(j);
                    result.(this_self+"_"+this_other)=d_fun(this_self,this_other);
                end
            end
            
            result=result{:,:}; %hack for now
        end
        
        
        function ges=runDETests(ges, tcounts, testname, options)
            arguments
                ges GroupedExpressionSummary
                tcounts {mustBeNumeric, mustBeEqualSizeCounts(tcounts,ges)}
                testname
                options
            end
            %perform a statistical DE test

            % multcompare correction
            % ctype='bonferroni'; %smooth values all the way to zero
            % % ctype='scheffe'; %smooth values all the way to zero

            % %genewise multiple comparison correction (for any test)
            % correctionmethod='fdr';
            
            switch testname
                case {'kw','kruskal','kruskalwallis'}
                case {'ranksum'}
                case {'anova1'}
                case {'ttest'}
                case {'chi2'} %proportions
                case {'ztest'} %proportions
                otherwise
                    error("unknown test: " + testname)
            end
        end
        
        
        function ges=findHVGs(ges, ncounts)
            % add HVG info
            
            % params.nBins
            % params.normMethod
            % params.minExpr
            % params.minCells
            % params.selectMethod
            % params.nHVG
            % params.dispThr
            %             
            %for plotting only:
%             genes.name=ges.name;
%             ,figID,highlightGenes
            
            [~, meanExpr, dispersion, rawnormdisp, ishvg]=findVariableGenes(ncounts,genes,params);
            
            ges.pooled.mean_n=meanExpr;
            ges.pooled.disp=dispersion;
            ges.pooled.norm_disp=rawnormdisp;
            ges.pooled.ishvg=ishvg;
        end
        
        function ges=subsetGenes(ges, subset)
            %remove genes by logical vector
        end
        
        function ges=subsetCells(ges, subset)
            %remove cells by logical vector - this changes the groups
            %allow new grouping to subset??
        end
        
    end

    methods (Access=private)
        function ges=computePairwiseContrasts(ges)
            
        end
    end
    
end

function mustBeEqualSizeCounts(a,ges)
    % Test for equal size
    if ~isequal(size(a),[ges.gene_count, ges.cell_count])
        eid = 'Size:notEqual';
        msg = 'Count matrix must be same size as the one used during construction';
        throwAsCaller(MException(eid,msg))
    end
end
        
function mustBeEqualRows(a,b)
    % Test for equal size
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notEqualRows';
        msg = 'Inputs must have equal number of rows.';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Inputs must have equal size.';
        throwAsCaller(MException(eid,msg))
    end
end

function groupLengthCheck(a,b)
    % Test for equal size
    if ~isequal(length(a),size(b,2))
        eid = 'Size:incorrectGroupLength';
        msg = 'group must be vector of length size(counts,2)';
        throwAsCaller(MException(eid,msg))
    end
end

function thrSizeCheck(a,b,c)
    % Test for correct # rows
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:incorrectThrRows';
        msg = 'thr must have same number of rows as counts';
        throwAsCaller(MException(eid,msg))
    end
    % Test for correct # cols
    if ~isequal(size(a,2),length(unique(c)))
        eid = 'Size:incorrectThrCols';
        msg = 'thr must have same number of columns as there are groups in thresh_group';
        throwAsCaller(MException(eid,msg))
    end
end

function thrVarNameCheck(a,b)
    %thr must be a table with varnames matching the thresh_group categories
    if ~isequal(sort(string(a.Properties.VariableNames(:))),sort(string(unique(b))))
        eid = 'Size:incorrectThrCols';
        msg = 'column names in thr must be the category names in thresh_group';
        throwAsCaller(MException(eid,msg))
    end
end

function groupNameCheck(a,b)
    %thr must be a table with varnames matching the thresh_group categories
    if ~ismember(string(a),string(b.group_names))
        eid = 'Name:incorrectGroupName';
        msg = 'column names in thr must be the category names in thresh_group';
        throwAsCaller(MException(eid,msg))
    end
end
