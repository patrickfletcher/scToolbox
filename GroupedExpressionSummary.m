classdef GroupedExpressionSummary
    %data class to compute and store gene expression data summarized for a
    %grouping variable

    %data is read only, must set values via functions (SetAccess=private)
    %need set access to allow sorting tables...
    properties 
        gene_count
        cell_count
        
        %pooled per-gene values (available without group) - one table?
        id %gene id
        name %gene symbol
        pooled = table() %n_cells, n_umi, stats, is_hvg, ...
        
        %group-wise descriptive statistics tables
        % [nGenes x nGroups], rownames=gene_name, varnames=groupnames
        % access pattern: ges.mean_t.group1("gene_name") or ges.mean_t(gene_names,group_names)
        
        %grouping variable
        group %grouping variable (value = missing supported)
        group_ix %logical index for non-missing cells in group, into original counts matrices
        group_count
        group_names
        mean_t = table() %mean log normalized counts
        mean_n = table() %mean normalized counts (log of this is not same as above)
%         mean_rank %normalized to 0,1? this is actual effect size measured by kw
        prct = table() %>0
        asin_prop = table() %arcsine transformed fraction expressing
        
        
        %thresholds - one column per thresh_group category
        thresh_group %will use only same cells as group
        thresh_group_count
        thresh_group_names
        thr = table()
        prct_thr = table() %gene-wise threshold using Otsu's method (tcounts)
        asin_prop_thr = table()%using thr
        
        % statistical tests
        
        %global tests among groups
        p_anova double %columns - kw: mean ranks, one: anova1, means
        p_chi2 double %chi2 test for proportions
        
        %pairwise contrast statistics - tables, varname=test name?
        p_mc = table() %result from multcompare following anova
        p_z = table() %pairwise z-tests for proportions
        
        
        %pairwise effect sizes - compute all pairwise + lookup, vs use
        %functions: compute for requested self/other groups?
        %functions are pretty fast
%         fc_n
%         fc_t
%         lfc %problem: several subtle variations possible
%         d_prct
%         d_prct_thr
%         cohen_h
%         cohen_h_thr
    end
   
    
    methods
        function [ges, tcounts] = GroupedExpressionSummary(genes, counts, ncounts, tcounts, group, thresh_group, thr)
            arguments
                genes table
                counts {mustBeNumeric, mustBeEqualRows(counts,genes)}
                ncounts {mustBeNumeric, mustBeEqualSize(ncounts,counts)}
                tcounts {mustBeNumeric, mustBeEqualSize(tcounts,counts)}
                group {groupLengthCheck(group,counts)}
                thresh_group {groupLengthCheck(thresh_group,counts)} = ones(1,size(counts,2))
                thr {thrSizeCheck(thr,counts,thresh_group)} = zeros(1,size(counts,1))
            end
            %constructor.
            % -genes is a table containing the basics: id, name, n_cells, n_umi
            % -ncounts, tcounts: normalized, log1p transformed counts [ng,nc]
            % -group (optional): grouping variable - length nc vector (can have missing values)
            % -threshold_group + thr: Otsu thresholds must be computed on full raw
            % dataset, they are for ambient RNA...
            
            %populate basics
            ges.gene_count=size(counts,1); %ng
            ges.cell_count=size(counts,2); %nc
            
            ges.id=genes.id;
            ges.name=genes.name;
            
            %process groups if supplied
            ges=ges.groupBy(counts, ncounts, tcounts, group);
                
            %process thresholds if supplied
            if ~isequal(thresh_group, ones(1,size(counts,2))) && ~isequal(thr, zeros(1,size(counts,1)))
                [ges, tcounts]=ges.addThresholds(tcounts, thresh_group, thr);
            end
            
            ges=ges.computePairwiseContrasts();
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
            
            ges.group_ix=~ismissing(group);
            ges.group=group;
            ges.group_names=categories(group);
            ges.group_count=length(ges.group_names);
            
            tic
            %get the pooled stats restricted to valid cells
            ges.pooled.name=ges.name;
            ges.pooled.Properties.RowNames=ges.id;
            ges.pooled.n_cells=sum(counts(:,ges.group_ix)>0,2);
            ges.pooled.n_umi=sum(counts(:,ges.group_ix),2);
%             ges.pooled.mean_n=mean(ncounts,2);
%             ges.pooled.var_n=var(ncounts,[],2);
%             ges.pooled.disp=(ges.pooled.var_n-ges.pooled.mean_n)./(ges.pooled.mean_n).^2;
            ges.pooled.mean_t=mean(tcounts(:,ges.group_ix),2);
%             ges.pooled.var_t=var(tcounts(:,ges.group_ix),[],2);
            toc
            
            tic
            ges.mean_n.id=ges.id;
            ges.mean_t.id=ges.id;
            ges.prct.id=ges.id;
            ges.asin_prop.id=ges.id;
            ges.mean_n.name=ges.name;
            ges.mean_t.name=ges.name;
            ges.prct.name=ges.name;
            ges.asin_prop.name=ges.name;
            for i=1:ges.group_count
                thisName=ges.group_names{i};
                thisGroup=group==thisName;
                ges.mean_n.(thisName)=mean(ncounts(:,thisGroup),2);
                ges.mean_t.(thisName)=mean(tcounts(:,thisGroup),2);
                frac=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup);
                ges.prct.(thisName)=frac*100;
                ges.asin_prop.(thisName)=asin(sqrt(frac));
            end
            ges.mean_n.Properties.RowNames=ges.id;
            ges.mean_t.Properties.RowNames=ges.id;
            ges.prct.Properties.RowNames=ges.id;
            ges.asin_prop.Properties.RowNames=ges.id;
            toc
        end
        
        function [ges, tcounts]=addThresholds(ges, tcounts, thresh_group, thr)
            arguments
                ges GroupedExpressionSummary
                tcounts {mustBeNumeric, mustBeEqualSizeCounts(tcounts,ges)} 
                thresh_group {groupLengthCheck(thresh_group,tcounts)}
                thr table {thrSizeCheck(thr,tcounts,thresh_group), thrVarNameCheck(thr, thresh_group)}
            end
            %adds gene-wise Otsu thresholds for proportion measures
            % The thresholds should be computed with ALL cells in raw data!
            tic
            
            if ~iscategorical(thresh_group)
                thresh_group=categorical(thresh_group);
            end
            thresh_group(~ges.group_ix)=missing();
            thresh_group=removecats(thresh_group);
            
            ges.thresh_group=thresh_group;
            ges.thresh_group_names=categories(thresh_group);
            ges.thresh_group_count=length(ges.thresh_group_names);
            
            ges.thr=thr;
            ges.thr.name=ges.name;
            ges.thr.id=ges.id;
            ges.thr.Properties.RowNames=ges.id;
            ges.thr=movevars(ges.thr,'name','Before',1);
            
            THR=zeros(size(tcounts));
            for i=1:ges.thresh_group_count
                thisThrGroup=ges.thresh_group==ges.thresh_group_names{i};
                THR(:,thisThrGroup)=repmat(ges.thr.(ges.thresh_group_names{i}),1,nnz(thisThrGroup));
            end
            tcounts=tcounts-THR;
            
            ges.prct_thr.id=ges.id;
            ges.asin_prop_thr.id=ges.id;
            ges.prct_thr.name=ges.name;
            ges.asin_prop_thr.name=ges.name;
            for i=1:ges.group_count
                thisName=ges.group_names{i};
                thisGroup=ges.group==thisName;
                frac=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup);
                ges.prct_thr.(thisName)=frac*100;
                ges.asin_prop_thr.(thisName)=asin(sqrt(frac));
            end
            ges.prct_thr.Properties.RowNames=ges.id;
            ges.asin_prop_thr.Properties.RowNames=ges.id;
            
            toc
        end
        
        %return tables?? or just matrices?
        function result = getFoldChanges(ges, selfnames, othernames, genes, whichExpr, transform_method)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                genes = true(ges.gene_count,1)
                whichExpr = 'norm'
                transform_method = ''
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            if isempty(genes), genes=true(ges.gene_count,1); end
            
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
            
            result=result{:,:}; %hack for now
        end
        
        
        function h = getCohenH(ges, selfnames, othernames, genes, useThr)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                genes = true(ges.gene_count,1)
                useThr = false
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            if isempty(genes), genes=true(ges.gene_count,1); end
            
            if useThr
                h_fun=@(self,other) 2*(ges.asin_prop_thr{genes,self}-ges.asin_prop_thr{genes,other});
            else
                h_fun=@(self,other) 2*(ges.asin_prop{genes,self}-ges.asin_prop{genes,other});
            end
                    
            h=table();
            for i=1:length(selfnames)
                this_self=selfnames(i);
                for j=1:length(othernames)
                    this_other=othernames(j);
                    h.(this_self+"_"+this_other)=h_fun(this_self,this_other);
                end
            end
            
            h=h{:,:}; %hack for now
        end
        
        
        function prct_diff = getPrctDiff(ges, selfnames, othernames, genes, useThr)
            arguments
                ges GroupedExpressionSummary
                selfnames {groupNameCheck(selfnames, ges)}
                othernames {groupNameCheck(othernames, ges)}
                genes = true(ges.gene_count,1)
                useThr = false
            end
            selfnames=string(selfnames); %force strings
            othernames=string(othernames);
            if isempty(genes), genes=true(ges.gene_count,1); end
            
            if useThr
                d_fun=@(self,other) ges.prct_thr{genes,self}-ges.prct_thr{genes,other};
            else
                d_fun=@(self,other) ges.prct{genes,self}-ges.prct{genes,other};
            end
            
            prct_diff=table();
            for i=1:length(selfnames)
                this_self=selfnames(i);
                for j=1:length(othernames)
                    this_other=othernames(j);
                    prct_diff.(this_self+"_"+this_other)=d_fun(this_self,this_other);
                end
            end
            
            prct_diff=prct_diff{:,:}; %hack for now
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
                case {'anova1'}
                case {'proportions'}
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
