function hdp=dotplot(genes, ncounts, tcounts, glist, ident, ctnames, options)
arguments
    genes
    ncounts
    tcounts
    glist
    ident
    ctnames=[]
    options.Ggroup=[]
    options.sortby='size'
    options.figid=[]
    options.figunits='inches'
    options.figpos=[]
    options.margins=[0.1,0.1,0.1,0.1]
    options.gap=0.1
    options.area_range=[1,144]
    options.min_dot_prct = 5
    options.min_any_prct = 5
    options.min_any_ident = []
    options.min_all_prct=0
    options.min_all_ident = []
    options.max_all_prct=100
    options.max_all_ident = []
    options.max_any_prct=100
    options.max_any_ident = []
    options.prct_leg=3
    options.prct_leg_height=0.05
    options.cb_prct_gap=0.01
    options.cb_gap=0.02
    options.cb_width=0.02
    options.cblabel='log mean expression'
    options.expr_method='mean'
    options.only_expressing=false
    options.verbose=false
    options.log=true
    options.logbase='e'
    options.pre_scale=false
    options.scale_method=[]
    options.scale_param=[]
    options.is_diverging=false
    options.cmap=[];
    options.sortgroups='none'
end

%TODO: refator plotDotPlot

if ~isempty(options.scale_method)
    if ismember(options.scale_method,["zscore","center","medianiqr"])
        options.is_diverging=1;
    end
    switch options.scale_method
        case 'zscore'
            options.cblabel="z score";
        case 'center'
            options.cblabel="centered expression";
        case {'norm','range'}
            options.cblabel="normalized expression";
    end
end

options.marg_h=options.margins(1:2);
options.marg_w=options.margins(3:4);
options.minArea=options.area_range(1);
options.maxArea=options.area_range(2);

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

ident=removecats(ident);
identnames=categories(ident);
cellsub=~ismissing(ident);
ident=ident(cellsub);
if isempty(ctnames)
    ctnames=identnames;
end

[gix, G]=getGeneIndices(glist,genes.name);

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

[exprTab,EXPR,PRCT]=getExpression(genes(gix,:),N,T,ident,method=options.expr_method,only_expressing=options.only_expressing);

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

%restrict to CTs to be plotted
[~,ctix]=ismember(ctnames,identnames);
EXPR=EXPR(:,ctix);
PRCT=PRCT(:,ctix);

discard=lowpct|not_all_min|all_expr|any_expr;
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

%hack until revamp plotDotPlot
options.min_prct=options.min_dot_prct;

if options.log
    switch options.logbase
        case '2'
            EXPR=log2(EXPR+1);
        case '10'
            EXPR=log10(EXPR+1);
        case 'e'
            EXPR=log1p(EXPR);
    end
end

if ~options.pre_scale
    if ~isempty(options.scale_method)
        if isempty(options.scale_param)
            EXPR=normalize(EXPR,2,options.scale_method);
        else
            EXPR=normalize(EXPR,2,options.scale_method, options.scale_param);
        end
    end
end

Glab=strcat('\it',G);
[ax,hs,cb,varNames]=plotDotPlot(Glab,ctnames,PRCT,EXPR,fh,Ggrp,options.sortby,options);

if isempty(options.cmap)
    if options.is_diverging
        cmap=split_cmap();
    else
        cmap=cbrewer('seq','Reds',64); cmap=[1,1,1;cmap];
    end
end
colormap(ax(1),cmap);

if options.is_diverging
    ax(1).CLim=max(abs(EXPR(:)))*[-1,1];
    cb.Limits=[min(EXPR(:)),max(EXPR(:))];
end

varNames=strrep(varNames,"\it","");
hdp.ax=ax;
hdp.hs=hs;
hdp.cb=cb;
exprTab.discard=discard(:);
[~,locB]=ismember(glist,varNames);
exprTab.rowperm=locB(:);
hdp.exprTab=exprTab;
hdp.options=options;
