classdef GroupedExpressionSummary < handle
    %data class to compute and store gene expression data summarized for a
    %grouping variable. Gene-wise summary quantities are computed across
    %all cells (pooled) and per group (grouped). Pairwise comparisons can
    %also be done, along with statistical tests, and combined per group to
    %give a summary.  Inspired by scran "scoreMarkers".

    %Group summaries
    % - basic summary of expression distribution per group
    % https://www.mathworks.com/help/matlab/descriptive-statistics.html
    % https://www.mathworks.com/help/stats/descriptive-statistics.html
    %
    % ideas: 
    % - bi/multimodality indicator: should this group be split?
    % - stats toolbox: grpstats

    % Support for blocking variable
    % - add "block" level in structures: usual analysis, but per block

    %Effect sizes
    % https://en.wikipedia.org/wiki/Effect_size 
    % 1. self vs others(pooled): use grand mean of other groups' means
    % - this allows each other group to contribute equally regardless of group sizes
    % 2. pairwise: evaluate every pairwise effect, summarize (scoreMarkers)
    % ideas: 
    % - differences in central tendency...
    % - correlation-type effect sizes...

    %tests between distributions
    % - statistics are also effect sizes (e.g., f=U/(n1*n2)=AUC in ranksum)
    % - also provide p-values (use FDR)
    % ideas:
    % - bootstrap nulls?

    %Filters
    % - Apply to pooled, grouped, and pairwise quantities
    % - combine filters to extract final gene sets
    % ideas:



    % in SCDataset, I can wrap this class with SCD methods
    %scd would already have: gene info, n/tcounts, cell metadata.  could
    % -> just pass in a grouping variable + options


    %TODO:
    % - handle nan/inf, small/missing group (or group-per-block)

    % - if we transpose expr, can use simpler functions
    % (@mean instead of @(x)mean(x,2).  Just transpose the
    % result back for table storage...

    %AUC: just use ranksum.  Similarly, could implement any of the relevant
    %2-sample tests if they give a stat that can be used as effect size. 

    % p-vals are not necessarily biologically meaningful. Would prefer some
    % kind of bootstrap methods: generate nulls from the data


    properties

        options

        %pooled per-gene values (available without group) - one table?
        id %gene id
        gene %gene symbol
        n_genes
        n_cells
        cellsub
        genesub

        %grouping variable
        group %grouping variable (value = missing supported)
        n_groups
        group_names %cols of per-group matrices
        group_counts

        %blocking variable.
        % https://rdrr.io/github/LTLA/scuttle/man/correctGroupSummary.html
        block
        n_blocks
        block_names
        block_counts

        %store sparse expression matrix? 
        % - slow ops when matrix is small. Make full within ops?
        % - remove on save, regenerate on load?? 
        % -- datafile + full_subset giving group (not bad)
        expr

        % pooled items across all cells. add arbitrary items easily to table
        pooled table %n_cells, n_umi, stats, is_hvg, ...

        %group-wise quantities (mean, prop, ...)
        % needs the batch correction using "block"
        grouped

        %pairwise effect sizes
        %store a struct, fields are cell types. within, a table for each
        %effect size can be added. Also add summary table like
        %"scoreMarkers": mean/median/min/max/rank per effectsize
        effectsizes

        % store logicals on pooled, grouped, or pairwise quantities
        % filters

%         %place to store parameters used?
%         expr_thr
%         min_cells
%         transform
%         invtransform
% 
%         % store params used on any additional function calls? by key?
%         params
    end


    methods
        % for simplicity, stick to one representation of expression at a
        % time. lognormcounts like scoreMarkers
        function ges = GroupedExpressionSummary(genes, expr, group, options)
            arguments
                genes table
                expr {mustBeNumeric, mustBeEqualRows(expr,genes)}
                group {groupLengthCheck(group,expr)}
                options.block=[]
                options.genesub=[]
                options.expr_thr=0
                options.expr_thr_method="prctile"
                options.min_cells=0
                options.block_reduce="grand_mean"
                options.transform="" %specify what transform was done
                options.invtransform="" %and how to reverse it if needed?
                options.only_expressing=false
                options.nanflag {mustBeMember(options.nanflag,["includenan","omitnan"])} ='omitnan'
                options.verbose=false
            end
            %constructor.
            % -genes is a table containing the basics: id, name, n_cells, n_umi
            % -expr: counts, ncounts, tcounts [ng,nc]
            % -group (optional): grouping variable - length nc vector (can have missing values)
            % -block: use blocking?
            % -- correct using linear model of group summary stats like scran?
            
            %group setup
            if ~iscategorical(group)
                if isnumeric(group)
                    group = compose("c%.2d",group);
                end
                group = categorical(group);
            end
            group = removecats(group(:));

            %block setup
            block=options.block;
            if isempty(block)
                block=ones(1,size(expr,2));
            end
            if ~iscategorical(block)
                block=categorical(block);
            end
            block=removecats(block(:));

            % cell subsetting to non-missing observations in group
            ges.cellsub = ~ismissing(group); 
            if ~isequal(nnz(ges.cellsub),size(expr,2))
                expr = expr(:,ges.cellsub);
                group = group(ges.cellsub);
                block = block(ges.cellsub);
            end

            %store group and block properties
            ges.group=group;
            ges.group_names=string(categories(group));
            ges.group_counts=countcats(group);
            ges.n_groups=length(ges.group_names);
            ges.block=block;
            ges.block_names=string(categories(block));
            ges.block_counts=countcats(block);
            ges.n_blocks=length(ges.block_names);
            
            %process genes: can remove genes not detected in this cell subset
            if isempty(options.genesub) 
                options.genesub = true(size(expr,1),1);
            end
            ges.genesub=options.genesub;
            cells_per_gene = sum(expr > 0,2);
            low_detected = cells_per_gene < options.min_cells;
            ges.genesub = ges.genesub & ~low_detected;

            if ~isequal(nnz(ges.genesub),size(expr,1))
                genes = genes(ges.genesub,:);
                expr = expr(ges.genesub,:);
            end

            ges.options=options;

            ges.id=genes.id;
            ges.gene=make_vars_unique(genes.name); %needs makeVarsUnique
            [ges.n_genes, ges.n_cells]=size(expr);

            %storing matrix?
%             ges.expr=sparse(expr);

            %get the pooled stats restricted to valid cells
            tic
            ges.pooled.id=ges.id;
            ges.pooled.name=ges.gene;
            ges.pooled.Properties.RowNames=ges.gene;
            ges.computeSummaries(expr, ["mean","prop","std"], use_group=0);
            
            if options.verbose
            toc
            end

            %compute group summaries
            if ges.n_groups>1
                tic
                % basic grouped values
                ges.computeSummaries(expr, ["mean","prop","std"]);
                
                if options.verbose
                toc
                end

                % initialize the effectsize struct
                for i=1:length(ges.group_names)
                    self=ges.group_names(i);
                    ges.effectsizes.(self).summary=table();
                    ges.effectsizes.(self).summary.gene=ges.gene;
                end
% %TODO: API - support passing list of multiple summary names to be computed
%                 ges.pairwiseDeltas("fc_expr");
%                 ges.pairwiseDeltas("cohen_d");
%                 ges.pairwiseDeltas("delta_prop");
%                 ges.pairwiseDeltas("lfc_prop");
%                 ges.pairwiseDeltas("cohen_h");
%                 toc
            end
        end

        
        %computeSummaries entry function
        % - support list of named methods
        % - support using group/block or not: 
        %   {all pooled, grouponly (blocks pooled), blockonly (groups pooled), both (groups, block-corrected)}
        % - if not specified, default keys used: method name(s)
        % - if specified, must provide one key per method
        % - method can be a function handle: method=@(x) ... key required 
        % - some common extra arguments can be passed as name/val pairs
        function result = computeSummaries(ges, X, methods, keys, options)
            arguments
                ges
                X
                methods = "all"
                keys = ""
                options.use_group=1
                options.use_block=1
                options.block_reduce="grand_mean"
            end

            [keys, funcs] = parseSummaryMethods(methods, keys, ges.options);

            if ges.options.only_expressing
                X(X==0)=nan;
            end

            groups=categorical(ones(size(ges.group)));
            if options.use_group
                groups=ges.group;
            end
            gnames=categories(groups);
            nG=length(gnames);

            blocks=categorical(ones(size(ges.group)));
            if options.use_block
                blocks=ges.block;
            end
            bnames=categories(blocks);
            nB=length(bnames);


            for ii=1:length(funcs)
                method_fcn=funcs{ii};

                result=zeros(size(X,1),nG*nB);
                for i=1:nG
                    for j=1:nB
                        thisGroup=groups==gnames(i) & blocks==bnames(j);
                        x=X(:,thisGroup);
                        result(:,(i-1)*nB+j) = method_fcn(x);
                    end
                end

                if nB>1
                    switch options.block_reduce
                        case "grand_mean"
                            result=ges.averageBlocks(result, groups, blocks,"mean");
                        case "linear"
                            result=ges.correctGroupSummary(result, groups, blocks);
                        case "none"
                            groups=stratify_factors(groups,blocks);
                            gnames=categories(groups);
                    end
                end

                if options.use_group
                    result=array2table(result, VariableNames=gnames, RowNames=ges.gene);
                    ges.grouped.(keys(ii))=result;
                else
                    ges.pooled.(keys(ii))=result;
                end
            end
        end

        function result = averageBlocks(ges, X, group, block, reduce_fcn)
            arguments
                ges
                X %group summaries
                group = ges.group
                block = ges.block
                reduce_fcn="mean"
            end
            nG=length(categories(group));
            nB=length(categories(block));

            switch reduce_fcn
                case "mean"
                    fcn=@(x)mean(x,2,"omitnan");
                case "median"
                    fcn=@(x)median(x,2,"omitnan");
            end
            result=zeros(size(X,1),nG);
            for i=1:nG
                result(:,i)=fcn(X(:,(i-1)*nB+(1:nB)));
            end
        end
        
        function result = correctGroupSummary(ges, X, group, block, options)
            arguments
                ges
                X %group summaries
                group = ges.group
                block = ges.block
                options.transform = "raw"
                options.offset = []
                options.weights = []
            end

            %transform the group summaries
            switch options.transform
                case "log"
                    if isempty(options.offset)
                        X = log1p(X);
                    else
                        X = log(X + options.offset);
                    end
                case "logit"
                    if isempty(options.offset)
                        options.offset=0.01;
                    end
                    X = (X + options.offset) / (1 + 2 * options.offset);
                    X = log(X./(1-X));
            end
            
            gnames=string(categories(group));
            nG=length(gnames);
            bnames=string(categories(block));
            nB=length(bnames);

            gind=categorical(repelem(gnames,nB,1),gnames,gnames);
            bind=categorical(repmat(bnames,nG,1),bnames,bnames);
            gtab_base=table();
%             gtab_base.group=gind;
%             gtab_base.block=bind;

            gtab_base=[gtab_base,array2table(dummyvar(gind),VariableNames=gnames)];
            db=dummyvar(bind);
            gtab_base=[gtab_base,array2table(db(:,2:end),VariableNames=bnames(2:end))];

            pnames=gtab_base.Properties.VariableNames;
% R:
% g=factor(c('g1','g1','g2','g2','g3','g3'))
% b=factor(c('b1','b2','b1','b2','b1','b2'))
% model.matrix(~0+g+b)
%   gg1 gg2 gg3 bb2
% 1   1   0   0   0
% 2   1   0   0   1
% 3   0   1   0   0
% 4   0   1   0   1
% 5   0   0   1   0
% 6   0   0   1   1

            result=zeros(size(X,1), nG);
            for i=1:size(X,1)
                thistab=gtab_base;
                thistab.y=X(i,:)';
%                 mdl=fitlm(thistab,"y~-1+block+group");

                %NOTES
                % - mdl here always seems to exclude the first group (Le,P1)
                % - use Coeff? residuals? as result

                mdl=fitlm(thistab,ResponseVar="y",PredictorVars=pnames,Intercept=false);
                Y(i,:)=mdl.Coefficients.Estimate(1:nG);
            end

            %reverse the transformation
            switch options.transform
                case "log"
                    if isempty(options.offset)
                        Y = expm1(Y);
                    else
                        Y = exp(Y) - options.offset;
                    end
                    Y(Y<0)=0;

                case "logit"
                    if isempty(options.offset)
                        options.offset=0.01;
                    end
                    Y = exp(Y);
                    Y = Y./(Y+1);
                    Y = Y*(1 + 2 * options.offset) - options.offset;
                    Y(Y<0)=0;
                    Y(Y>1)=1;
            end
            result=Y;
        end


        % Effect sizes
        % - some operate just on group summary quantities (already computed)
        % - others require original expression values:
        % -- pooled std of pair (cohen_d)
        % -- ranking in pooled values of pair (AUC)
        %=> optionally accept expression matrix as input?

        % compute effect sizes: compose as difference(transform(X))
        function [summary, all_stats] = pairwiseDeltas(ges, effect_name, selfnames, othernames, key, options)
            arguments
                ges GroupedExpressionSummary
                effect_name 
                selfnames = ges.group_names
                othernames = ges.group_names
                key = effect_name
                options.pre_transform function_handle = @(x)x
                options.post_transform function_handle = @(x)x
                options.offset=0.01
                options.save_pw_stats=true
                options.overwrite_key=false
            end
            % use for lfc in mean expression or proportions. e.g.,
            % - lfc_mean: log_mean(self)-log_mean(other)
            %             log(mean(self))-log(mean(other))  [using pretransform=log on mean normcounts]
            % - cohen_d (weighting s equally as done in score markers):
            %      [log_mean(self)-log_mean(other)]/mean(std(self),std(other))
            %
            % - lfc_prop: log(prop(self)) - log(prop(other))
            % - delta_prop: prop(self) - prop(other)
            % - cohen_h: 0.5(asin(sqrt(prop(self))) - 0.5(asin(sqrt(prop(other)))
            % - delta_prct: =100*delta_prop

            %could use self/othernames+key to just do specific contrasts

            %TODO: what if run same effectsize twice?
%             if isfield(ges.effectsizes, key) && options.overwrite_key==false
%                 warning("Key exists:" + key + ". Set overwrite_key=true to force overwrite.")
%                 return
%             end

            pre_transform=options.pre_transform;
            post_transform=options.post_transform;

            selfnames=string(selfnames); %force strings
            othernames=string(othernames);

            es_fun=@(x1, x2) x1-x2;
            switch effect_name
                case {"delta_expr","lfc_expr","cohen_d"}
                    %set options.pre_transform=@log1p if using normcounts
                    %will manually standardize for cohen_d
                    X=ges.grouped.mean{:,:};
                    post_transform = @(x)x; %force no post
                case "fc_expr"
                    X=ges.grouped.mean{:,:};
                    pre_transform = @(x)x+options.offset;
                    es_fun=@(x1,x2) x1./x2;
                case "delta_prop"
                    X=ges.grouped.prop{:,:};
                    post_transform = @(x)x; %force no post
                case "lfc_prop"
                    X=ges.grouped.prop{:,:};
                    pre_transform = @log1p;
                    post_transform = @(x)x; %force no post
                case "cohen_h"
                    X=ges.grouped.prop{:,:};
                    pre_transform = @(x)asin(sqrt(x));
                    post_transform = @(x)0.5*x;
                otherwise
                    disp("unknown/not implemented:" + effect_name)
            end

            X=pre_transform(X);

            for i=1:length(selfnames)
                self=selfnames(i);
                all_stats=ges.effectsizes.(self);
                selfix=ges.group_names==self;
                others=setdiff(othernames,self);
                otherix=ismember(ges.group_names,others);
                this_fc= es_fun(X(:,selfix),X(:,otherix));
                
                if effect_name=="cohen_d"
                    S=ges.grouped.std{:,otherix};
                    %broadcast self-std to size of others
                    S(:,:,2)=repmat(ges.grouped.std{:,selfix},1,length(others));
                    this_fc=this_fc./mean(S,3);
                end

                this_fc=post_transform(this_fc);     

%                 fc_tab=table();
%                 fc_tab.gene=ges.gene;
%                 fc_tab=[fc_tab,array2table(this_fc,"VariableNames",others)];
                fc_tab=array2table(this_fc,"VariableNames",others, "RowNames",ges.gene);
                all_stats.(key)=fc_tab;
                summary=ges.summarizeContrasts(fc_tab, key);
                all_stats.summary=[all_stats.summary,summary];
                ges.effectsizes.(self)=all_stats;
            end

        end

        % Aggregation across pairwise comparisons
        % - min/median/max, mean, minrank, ...
        function summary=summarizeContrasts(ges, this_ES, key, options)
            arguments
                ges
                this_ES
                key
                options.rankdir='descend'
            end

            summary=table();
            summary.("mean_"+key)=mean(this_ES{:,:},2);
            summary.("median_"+key)=median(this_ES{:,:},2);
            summary.("min_"+key)=min(this_ES{:,:},[],2);
            summary.("max_"+key)=max(this_ES{:,:},[],2);

            esRank=zeros(size(this_ES));
            for i=1:size(this_ES,2)
                [~,ixs]=sort(this_ES{:,i},options.rankdir,'MissingPlacement','last');
                esRank(ixs,i)=1:size(this_ES,1);
            end
            summary.("minrank_"+key)=min(esRank,[],2);
        end

        %wrapper function to support lists of effect sizes. If effect_names
        %not specified, compute a standard set.
        function computeEffectSizes(ges, effect_names, selfnames, othernames, options)
            arguments
                ges GroupedExpressionSummary
                effect_names = ["lfc_expr","cohen_d","delta_prop"]
                selfnames = ges.group_names
                othernames = ges.group_names
                options.save_pw_stats=true
                options.overwrite_key=false
            end

%             pwargs=namedargs2cell(options);

            for i=1:length(effect_names)
                pairwiseDeltas(ges, effect_names(i), selfnames, othernames);
            end
        end

        % extract the top markers across all clusters based on thresholding
        % effect size summary tables
        % TODO - do this simultaneously for multiple summary features? 
        function result = topMarkers(ges, summary_name, method, selfnames, othernames, options)
            arguments
                ges
                summary_name = "minrank_delta_prop"
                method = "thr"
                selfnames = ges.group_names
                othernames = ges.group_names
                options.nTop=10;
                options.thr=5;
                options.p=95;
                options.min_self_prop=0
                options.max_other_prop=1
                options.summary_min=-inf;
                options.summary_max=inf;
            end
            %default: keep minrank_delta_prop<=5

            sortdir="descend";
            thr_dir="gt";
            if contains(summary_name,"minrank")
                sortdir="ascend";
                thr_dir="lt";
            end

            result=table;
            for i=1:length(selfnames)
                thisname=selfnames(i);
                thisothers=setdiff(othernames,thisname);
                thistab=ges.effectsizes.(thisname).summary;

                minpropfilt=ges.grouped.prop.(thisname)>=options.min_self_prop;
                maxpropfilt=all(ges.grouped.prop{:,thisothers}<=options.max_other_prop,2);
                propfilt=minpropfilt&maxpropfilt;
                if nnz(propfilt)>0
                    thistab=thistab(propfilt,:);
                else
                    disp("No top markers satisfying proportion constraints: "+ thisname)
                    continue
                end
    
                thistab=sortrows(thistab,summary_name,sortdir,'missingplacement','last');
                keep = true(size(thistab,1),1);
                switch method
                    case "thr"
                        keep=thistab.(summary_name)>=options.thr; 
                        if thr_dir=="lt"
                            keep=thistab.(summary_name)<=options.thr; 
                        end
                    case "top"
                        keep = ismember(1:height(thistab), 1:min(options.nTop,height(thistab)))';
                    case "prctile"
                        thr=prctile(thistab.(summary_name),options.p);
                        keep=thistab.(summary_name)>=thr; 
                        if thr_dir=="lt"
                            keep=thistab.(summary_name)<=thr; 
                        end
                end

                %apply summary min/max
                keep = keep & thistab.(summary_name)>=options.summary_min;
                keep = keep & thistab.(summary_name)<=options.summary_max;

                if nnz(keep)>0
                    thistab=thistab(keep,:);
                    thistab.celltype=repmat(thisname,nnz(keep),1);
                    thistab=movevars(thistab,"celltype","before",1);
                    result=[result;thistab];
                else
                    disp("No top markers satisfying effect size constraints: "+ thisname)
                end
            end

        end

        % stats from hypothesis test functions
        % generate both effect sizes and p-values...
        % - ttest2 VarType=unequal (Cohen's d with pooled var?)
        % - ranksum (AUC)
        % - kstest2??
        function pairwiseTest(ges, testname)
            arguments
                ges
                testname
            end
        end


        %combine p-vals from statistical tests..
        function summary=combinePvals(ges, key, options)
            arguments
                ges
                key
                options
            end

        end


        % Select markers based on thresholds of effect sizes
        function result = findMarkers(ges,selfnames,othernames,options)
            arguments
                ges
                selfnames = ges.group_names
                othernames = ges.group_names
                options.effect_names="fc_expr"
                options.thresholds=2
                options.limits=1
                options.directions="up"
                options.effect_logic="and" %@(x)x %eg: @(x) (x(:,1)|x(:,2))&x(:,3)
                options.self_combine='all'
                options.self_minprop=0.25
                options.self_maxprop=1
                options.other_combine='all'
                options.other_minprop=0
                options.other_maxprop=1
            end

            selfnames=string(selfnames);
            othernames=string(othernames);

            nS=length(selfnames);
            nES=length(options.effect_names);
            es=options.effect_names;
            thr=options.thresholds;
            dir=options.directions;
            effect_logic=options.effect_logic;

            for i=1:nS
                self=selfnames(i);
                selfix=ges.group_names==self;
                others=setdiff(othernames,self);
                nO=length(others);
                otherix=ismember(ges.group_names,others);

                %self-group filters
                self_filt=ges.grouped.prop{:,selfix}>=options.self_minprop & ges.grouped.prop{:,selfix}<=options.self_maxprop;

                %other-group filter
                other_filt=ges.grouped.prop{:,otherix}>=options.other_minprop & ges.grouped.prop{:,otherix}<=options.other_maxprop;

                %pairwise ES filters
                es_filters=true(ges.n_genes,nES,nO);
                for j=1:nES
                    this_es=ges.effectsizes.(self).(es(j));
                    switch dir(j)
                        case "up"
                            this_es_filt=this_es{:,otherix}>=thr(j);
                        case "down"
                            this_es_filt=this_es{:,otherix}<=thr(j);
                    end

                    %broadcast other_filt to each pairwise es filter
                    this_es_filt=this_es_filt & other_filt;

                    es_filters(:,j,:)=this_es_filt;
                end

                %combine es filters into one
%                 combined_es = effect_logic(combined_pairwise);
                switch effect_logic
                    case "and"
                        combined_pairwise=all(es_filters,2);
                    case "or"
                        combined_pairwise=any(es_filters,2);
                end

                %combine the filters pairwise
                switch options.other_combine
                    case "all"
                        combined_es=all(combined_pairwise,3);
                    case "any"
                        combined_es=any(combined_pairwise,3);
                end

                
                filt=self_filt&combined_es;

                tab=table;
                result.(self).markers=ges.gene(filt);
                result.(self).effectsizes=ges.effectsizes.(self).summary(filt,:);
                result.(self).mean=ges.grouped.mean(filt,:);
                result.(self).prop=ges.grouped.prop(filt,:);

            end
        end


        % may need some util functions to extract/combine contrasts


    end

    methods (Access=private)
    end

end

function mustBeEqualSizeExpr(a,ges)
% Test for equal size
if ~isequal(size(a(ges.genesub,ges.cellsub)),[ges.n_genes, ges.n_cells])
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


function groupNameCheck(a,b)
%thr must be a table with varnames matching the thresh_group categories
if ~ismember(string(a),string(b.group_names))
    eid = 'Name:incorrectGroupName';
    msg = 'column names in thr must be the category names in thresh_group';
    throwAsCaller(MException(eid,msg))
end
end


%similar to "groupsummary"
function [keys, funcs] = parseSummaryMethods(methods, keys, options)
arguments
    methods
    keys
    options
end
    allMethods = {'mean', 'sum', 'min', 'max', 'range', 'median','var', 'std'};
    if nargin==0 && nargout<2
        keys=string(allMethods);
        return
    end

    switch options.expr_thr_method
        case "prctile"
            thr_fun = @(x) x > prctile(x,options.expr_thr);
        otherwise
            thr_fun = @(x) x > options.expr_thr;
    end

    %methods could be a cell array with mix of names/fcn handles. force
    %cellstr if it is a string array
    if isstring(methods)
        methods = cellstr(methods);
    elseif ~iscell(methods)
        methods = {methods};
    end
    numMethods = numel(methods);

    if methods=="all"
        methods=allMethods;
    end

    if length(keys)~=length(methods)
        keys = strings(1,numMethods);
    end

    numfun = 1;
    funcs = cell(1,numMethods);
    removeix=false(1,numMethods);
    for i = 1:numMethods
        if ischar(methods{i})
            switch methods{i}
                case 'mean'
                    funcs{i}=@(x) mean(x,2,options.nanflag);
                case 'median'
                    funcs{i}=@(x) median(x,2,options.nanflag);
                case 'min'
                    funcs{i}=@(x) min(x,[],2);
                case 'max'
                    funcs{i}=@(x) max(x,[],2);
                case 'ptile'
                    funcs{i}=@(x) prctile(x,options.ptile,2);
                case 'trimean'
                    funcs{i}=@(x) mean(prctile(x,[25,50,50,75],2),2);
                case 'trimmean'
                    funcs{i}=@(x) trimmean(x, options.trim_prct, 2);
                case 'iqr'
                    funcs{i}=@(x) iqr(x,2);
                case 'std'
                    funcs{i}=@(x) std(x,0,2);
                case 'std1'
                    funcs{i}=@(x) std(x,1,2);
                case 'var'
                    funcs{i}=@(x) var(x,0,2);
                case 'var1'
                    funcs{i}=@(x) var(x,1,2);
                case 'cv'
                    funcs{i}=@(x) std(x,0,2)./mean(x,2);
                case 'disp'
                    funcs{i}=@(x) (var(x,0,2)-mean(x,2))./mean(x,2).^2;
                case 'sum'
                    funcs{i}=@(x) sum(x,2);
                case 'num'
                    funcs{i}=@(x) sum(thr_fun(x),2);
                case 'prop'
                    funcs{i}=@(x) sum(thr_fun(x),2)./size(x,2);
                case 'prct'
                    funcs{i}=@(x) sum(thr_fun(x),2)./size(x,2)*100;
                otherwise
                    disp("method '"+methods{i}+"' not implemented")
                    removeix(i)=1;

            end
            if keys(i)==""
                keys(i)=string(methods{i});
            end
            
        else
            if ~isa(methods{i},'function_handle')
                error(message('MATLAB:groupsummary:InvalidMethodOption'));
            end
            if keys{i}==""
                keys = "fun"+num2str(numfun);
                numfun = numfun +1;
            end
        end
    end
    keys(removeix)=[];
    funcs(removeix)=[];
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function thrSizeCheck(a,b,c)
% % Test for correct # rows
% if ~isequal(size(a,1),size(b,1))
%     eid = 'Size:incorrectThrRows';
%     msg = 'thr must have same number of rows as counts';
%     throwAsCaller(MException(eid,msg))
% end
% 
% % Test for correct # cols
% if ~isequal(size(a,2),length(unique(c)))
%     eid = 'Size:incorrectThrCols';
%     msg = 'thr must have same number of columns as there are groups in thresh_group';
%     throwAsCaller(MException(eid,msg))
% end
% end

% function thrVarNameCheck(a,b)
% %thr must be a table with varnames matching the thresh_group categories
% if ~isequal(sort(string(a.Properties.VariableNames(:))),sort(string(unique(b))))
%     eid = 'Size:incorrectThrCols';
%     msg = 'column names in thr must be the category names in thresh_group';
%     throwAsCaller(MException(eid,msg))
% end
% end



%
% %special function for percent expressing
%         %                 expr %{mustBeEqualSizeExpr(expr,ges)}
%         function result = summarizeProportions(ges, X, thr, key, options)
%             arguments
%                 ges
%                 X
%                 thr = 0
%                 key = "prop"
%                 options.do_full=1
%             end
% %             X=ges.expr;
% %             if options.do_full
% %                 X=full(X);
% %             end
%             
%             groups=ges.group;
%             gnames=ges.group_names;
%             nG=ges.group;
%             gcounts=ges.group_counts;
%             if ges.n_blocks>1
%                 % first get summaries for groups x blocks
%                 groups=stratify_factors(ges.group,ges.block); %TODO: remove dependency
%                 gnames=categories(groups);
%                 gcounts=countcats(groups);
%                 nG=length(gcounts);
%             end
% 
%             result=zeros(ges.n_genes,nG);
%             for i=1:nG
%                 thisGroup=groups==gnames(i);
%                 x=X(:,thisGroup);
%                 result(:,i) = sum(x>thr,2)/gcounts(i);
%             end
% 
%             if ges.n_blocks>1
%                 result=ges.correctGroupSummary(result, transform="logit", offset=0.01);
%             end
% 
%             result=array2table(result,VariableNames=ges.group_names,RowNames=ges.gene);
% 
%             ges.grouped.(key)=result;
%         end
% 
% 
%         %                 expr %{mustBeEqualSizeExpr(expr,ges)}
%         function result = summarizeExpression(ges, X, method_fcn, key, options)
%             arguments
%                 ges
%                 X
%                 method_fcn function_handle
%                 key
%                 options.only_expressing = false
%                 options.pre_transform function_handle = @(x)x
%                 options.post_transform function_handle = @(x)x
%                 options.do_full = true
%             end
%             % populate the entries of the "grouped" struct. Provide an
%             % expression matrix, function handle applied to that matrix,
%             % and key to store it under.
%             %
%             % value associated with key table of values.
%             % pre/post transforms allow some flexibility...
%             % - use function handles only?
% 
% %             X=ges.expr;
% %             if options.do_full
% %                 X=full(X);
% %             end
% 
%             pre_transform=options.pre_transform;
%             post_transform=options.post_transform;
% 
% %TODO: this needs "omitnan" flag on common functions.
%             if options.only_expressing
%                 X(X==0)=nan;
%             end
% 
%             X=pre_transform(X);
% 
%             groups=ges.group;
%             gnames=ges.group_names;
%             nG=ges.group;
%             gcounts=ges.group_counts;
%             if ges.n_blocks>1
%                 % first get summaries for groups x blocks
%                 groups=stratify_factors(ges.group,ges.block); %TODO: remove dependency
%                 gnames=categories(groups);
%                 gcounts=countcats(groups);
%                 nG=length(gcounts);
%             end
% 
%             result=zeros(ges.n_genes,nG);
%             for i=1:nG
%                 thisGroup=groups==gnames(i);
%                 x=X(:,thisGroup);
%                 result(:,i) = method_fcn(x);
%             end
% 
%             %if specified, apply a post-transformation: e.g., log of mean(norm)
%             result=post_transform(result);
% 
%             if ges.n_blocks>1
%                 result=ges.correctGroupSummary(result);
%             end
% 
%             result=array2table(result,VariableNames=ges.group_names,RowNames=ges.gene);
%             ges.grouped.(key)=result;
%         end



%    function X = do_transform(X, method, options)
% arguments
%     X
%     method = ''
%     options.offset=1
% end
% %this is really for some of the special ones... listing others
% %for API consistency...
% switch method
%     case 'add_offset'
%         X = X + options.offset;
%     case 'log'
%         X = log(X);
%     case 'exp'
%         X = exp(X);
%     case 'log1p'
%         X = log1p(X);
%     case 'expm1'
%         X = expm1(X);
%         X(X<0)=0;
%     case 'log2p'
%         X = log2(X + options.offset);
%     case 'exp2m'
%         X = 2.^(X)-options.offset;
%         X(X<0)=0;
%     case 'log10p'
%         X = log10(X + options.offset);
%     case 'exp10m'
%         X = 10.^(X)-options.offset;
%         X(X<0)=0;
%     case 'log2'
%         X = log2(X);
%     case 'exp2'
%         X = 2.^(X);
%     case 'log10'
%         X = log10(X);
%     case 'exp10'
%         X = 10.^(X);
%     case 'logit'
%         X = (X + options.offset) / (1 + 2 * options.offset);
%         X = log(X./(1-X));
%     case 'rev_logit'
%         X = exp(X);
%         X = X./(X+1);
%         X = X*(1 + 2 * options.offset) - options.offset;
%         X(X<0)=0;
%         X(X>1)=1;
%     otherwise
%         %no op
% end
% end         
            
%             switch method
%                 case 'trimean'
%                     expr_fun=@(x) mean(prctile(x,[25,50,50,75],2),2);
%                 case 'trimmean'
%                     expr_fun=@(x) trimmean(x, options.trim_prct, 2);
%                 case 'sum'
%                     expr_fun=@(x) sum(x,2);
%                 case 'std'
%                     expr_fun=@(x) std(x,0,2);
% %                 case 'std1'
% %                     expr_fun=@(x) std(x,1,2);
%                 case 'var'
%                     expr_fun=@(x) var(x,0,2);
% %                 case 'var1'
% %                     expr_fun=@(x) var(x,1,2);
% %                 case 'cv'
% %                     expr_fun=@(x) std(x,0,2)./mean(x,2);
% %                 case 'disp'
% %                     expr_fun=@(x) (var(x,0,2)-mean(x,2))./mean(x,2).^2;
%                 otherwise 
%                     %catchall for common functions. safer to list them..
%                     expr_fun=eval("@(x)"+method+"(x,2,'omitnan')");
%             end



%         %generic getEffectSize? just do the basic ones...
%         function result = pairwiseEffectSize(ges, selfnames, othernames, effecttype,  slot, options)
%             arguments
%                 ges GroupedExpressionSummary
%                 selfnames {groupNameCheck(selfnames, ges)}
%                 othernames = ""
%                 effecttype = "lfc"
%                 slot = "mean_lognorm"
%                 options.genes = ""
%                 options.pre_transform=""
%                 options.post_transform=""
%                 options.standardize = true
%             end
%             selfnames=string(selfnames); %force strings
% 
%             if isempty(othernames)
%                 othernames = setdiff(ges.group_names,selfnames);
%             end
%             othernames=string(othernames);
% 
%             switch effecttype
%                 case "lfc"
%                     X=ges.grouped.(slot);
%                 case "AUC"
%                     %does this need the whole counts?
%                 case ["delta_prop","lfc_prop"]
% 
%             end
% %             X=ges.grouped.(slot);
% 
% 
% %             X=do_transform(X{gix,:}, options.pre_transform);
% %             if slot=="mean_norm" && ~contains(options.pre_transform,"log")
% %                 fc_fun=@(self,other) X(:,ges.group_names==self)./X(:,ges.group_names==other);
% %             else
% %                 fc_fun=@(self,other) X(:,ges.group_names==self)-X(:,ges.group_names==other);
% %             end
% 
% 
%             P=ges.grouped.prop;
%             d_fun=@(self,other) P{options.genes,self}-P{options.genes,other};
% 
%             for i=1:length(selfnames)
%                 this_self=selfnames(i);
%                 fc=d_fun(this_self, othernames);
%                 result.(this_self).(key)=array2table(fc,"RowNames",ges.gene(gix),"VariableNames",othernames);
%             end
%         end
% 
% 
%         %         function ges = correctGroupSummary(ges)
%         %         end



        %         function ges=runDETests(ges, tcounts, testname, options)
        %             arguments
        %                 ges GroupedExpressionSummary
        %                 tcounts {mustBeNumeric, mustBeEqualSizeCounts(tcounts,ges)}
        %                 testname
        %                 options
        %             end
        %             %perform a statistical DE test
        %
        %             % multcompare correction
        %             % ctype='bonferroni'; %smooth values all the way to zero
        %             % % ctype='scheffe'; %smooth values all the way to zero
        %
        %             % %genewise multiple comparison correction (for any test)
        %             % correctionmethod='fdr';
        %
        %             switch testname
        %                 case {'kw','kruskal','kruskalwallis'}
        %                 case {'ranksum'}
        %                 case {'anova1'}
        %                 case {'ttest'}
        %                 case {'chi2'} %proportions
        %                 case {'ztest'} %proportions
        %                 otherwise
        %                     error("unknown test: " + testname)
        %             end
        %         end


        %         function markers = findMarkers(ges, selfnames, othernames, options)
        %             arguments
        %                 ges GroupedExpressionSummary
        %                 selfnames {groupNameCheck(selfnames, ges)}
        %                 othernames {groupNameCheck(othernames, ges)}
        %                 options.direction = 'up'
        %                 options.selfpool = 'all'
        %                 options.otherpool = 'all'
        %                 options.fc_options = []  %pass a struct if desired
        %                 options.d_prct_options = []
        %                 options.fc_thr = 1
        %                 options.prct_diff = 0
        %                 options.min_self_prct = 20
        %                 options.max_self_prct = 100
        %                 options.min_other_prct = 0
        %                 options.max_other_prct = 5
        %                 options.anova_thr = 1
        %                 options.p_thr = 1
        %                 options.chi2_thr = 1
        %                 options.z_thr = 1
        %             end
        %             %apply conditions groupwise, pairwise to get list of genes
        %
        %             %self, other %
        %             self_prct=ges.prct(:,selfnames);
        %             other_prct=ges.prct(:,othernames);
        %
        %
        %             %get the effect sizes
        %             if isempty(options.fc_options)
        %                 fc=ges.getFoldChanges(selfnames,othernames);
        %             else
        %                 fc=ges.getFoldChanges(selfnames,othernames, options.fc_options);
        %             end
        %             if isempty(options.d_prct_options)
        %                 d_prct=ges.getPrctDiff(selfnames,othernames);
        %             else
        %                 d_prct=ges.getPrctDiff(selfnames,othernames, options.d_prct_options);
        %             end
        %
        %             %effect size selection
        %
        %
        %             %get DE p-values
        %             thisCombs=all(ismember(ges.pw_combs,selfnames)|ismember(ges.pw_combs,othernames),2); %index into columns of the P-matrices
        %
        %         end


        %         function ges=findHVGs(ges, ncounts)
        %             % add HVG info
        %
        %             % params.nBins
        %             % params.normMethod
        %             % params.minExpr
        %             % params.minCells
        %             % params.selectMethod
        %             % params.nHVG
        %             % params.dispThr
        %             %
        %             %for plotting only:
        % %             genes.name=ges.gene;
        % %             ,figID,highlightGenes
        %
        %             [~, meanExpr, dispersion, rawnormdisp, ishvg]=findVariableGenes(ncounts,genes,params);
        %
        %             ges.pooled.mean_n=meanExpr;
        %             ges.pooled.disp=dispersion;
        %             ges.pooled.norm_disp=rawnormdisp;
        %             ges.pooled.ishvg=ishvg;
        %         end