function hdp=dotplot(genes, ncounts, tcounts, glist, ident, ctnames, options, dotopts)
arguments
    genes
    ncounts
    tcounts
    glist
    ident
    ctnames=[]
    options.Ggroup=[]
    options.figid=[]
    options.figunits='inches'
    options.figpos=[]
    options.min_any_prct = 5
    options.min_any_ident = []
    options.min_all_prct=0
    options.min_all_ident = []
    options.max_all_prct=100
    options.max_all_ident = []
    options.max_any_prct=100
    options.max_any_ident = []
    options.max_sum_prct=100
    options.max_sum_ident = []
    options.expr_method='mean'
    options.expr_vals='ncounts'
    options.only_expressing=false
    options.verbose=false
    options.log=true
    options.logbase='10'
    options.pre_scale=false
    options.scale_method=[]
    options.scale_param=[]
    dotopts.sortby='size'
    dotopts.is_diverging=false
    dotopts.cmap=[];
    dotopts.min_dot_prct = 5
    dotopts.area_range=[1,144]
    dotopts.margins=[0.1,0.1,0.1,0.1]
    dotopts.gap=0.1
    dotopts.prct_leg=3
    dotopts.prct_leg_dims=[0.25,0.05]
    dotopts.cb_prct_gap=0.01
    dotopts.cb_gap=0.01
    dotopts.cb_dims=[0.3,0.01]
    dotopts.cb_digits=2;
    dotopts.sortgroups='none'
end

%TODO: refator plotDotPlot

switch options.expr_vals
    case 'ncounts'
        dotopts.cblabel='log mean expression';
    case 'tcounts'
        dotopts.cblabel='mean log expression';
end


if ~isempty(options.scale_method)
    if ismember(options.scale_method,["zscore","center","medianiqr"])
        options.is_diverging=1;
    end
    switch options.scale_method
        case 'zscore'
            dotopts.cblabel="z score";
        case 'center'
            dotopts.cblabel="centered expression";
        case {'norm','range'}
            dotopts.cblabel="normalized expression";
    end
end

if isempty(options.figid)
    fh=figure();
else
    fh=figure(options.figid);
end
fh.Units=options.figunits;
if ~isempty(options.figpos)
    fh.Position=options.figpos;
end

Ggrp=ones(size(glist));
if ~isempty(options.Ggroup) && length(options.Ggroup)==length(glist)
    Ggrp=options.Ggroup;
end
if ~iscategorical(ident)
    ident=categorical(ident);
end
ident=removecats(ident);
cellsub=~ismissing(ident);

ident=ident(cellsub);
ident=removecats(ident);
identnames=categories(ident);
if isempty(ctnames)
    ctnames=identnames;
end

[gix, G]=getGeneIndices(glist,genes.name);
Ggrp=Ggrp(ismember(glist,G));

N=ncounts(gix,cellsub);
T=tcounts(gix,cellsub);

if options.pre_scale
    if ~isempty(options.scale_method)
        if isempty(options.scale_param)
            N=normalize(N,2,options.scale_method);
        else
            N=normalize(N,2,options.scale_method, options.scale_param);
        end
    end
end

[exprTab,EXPR,PRCT]=getExpression(genes(gix,:),N,T,ident, ...
    expr_vals=options.expr_vals,method=options.expr_method, only_expressing=options.only_expressing);

%minimum prct in any cts of interest (ctnames by default)
if isempty(options.min_any_ident)
    options.min_any_ident=ctnames;
end
[~,minctix]=ismember(options.min_any_ident,identnames);
lowpct=all(PRCT(:,minctix)<options.min_any_prct,2);

%keep only those with min_all_prct in cts of interest (ctnames by default)
if isempty(options.min_all_ident)
    options.min_all_ident=ctnames;
end
[~,allctix]=ismember(options.min_all_ident,identnames);
not_all_min=any(PRCT(:,allctix)<options.min_all_prct,2);

%prct considered expressed for all-expressing exclusion in cts of interest (ctnames by default)
if isempty(options.max_all_ident)
    options.max_all_ident=ctnames;
end
[~,allctix]=ismember(options.max_all_ident,identnames);
all_expr=all(PRCT(:,allctix)>options.max_all_prct,2);

if isempty(options.max_any_ident)
    options.max_any_ident=ctnames;
end
[~,allctix]=ismember(options.max_any_ident,identnames);
any_expr=any(PRCT(:,allctix)>options.max_any_prct,2);

sum_prct=false(size(PRCT,1),1);
if ~isempty(options.max_sum_ident)
    [~,allctix]=ismember(options.max_sum_ident,identnames);
    sum_prct=sum(PRCT(:,allctix),2)>options.max_sum_prct;
end

%restrict to CTs to be plotted
[~,ctix]=ismember(ctnames,identnames);
EXPR=EXPR(:,ctix);
PRCT=PRCT(:,ctix);

discard=lowpct|not_all_min|all_expr|any_expr|sum_prct;
EXPR(discard,:)=[];
PRCT(discard,:)=[];
if options.verbose && nnz(discard)>0
    disp("Genes expressed in less than min_prct for all cell groups:")
    disp(G(lowpct)')
    disp("Genes expressed in more than all_expr_prct for all cell groups:")
    disp(G(all_expr)')
end
G(discard)=[];
Ggrp(discard)=[];

if options.expr_vals=="ncounts" && options.log
    switch options.logbase
        case '2'
            EXPR=log2(EXPR+1);
        case '10'
            EXPR=log10(EXPR+1);
        case 'e'
            EXPR=log1p(EXPR);
    end
end

if options.pre_scale
    if ~isempty(options.scale_method)
        if isempty(options.scale_param)
            EXPR=normalize(EXPR,2,options.scale_method);
        else
            EXPR=normalize(EXPR,2,options.scale_method, options.scale_param);
        end
    end
end

Glab=strcat('\it',G);

%hack until revamp plotDotPlot
dotargs=namedargs2cell(dotopts);
hdp=plotDotPlot(Glab,ctnames,PRCT,EXPR,fh,Ggrp,dotargs{:});

hdp.options=options;
