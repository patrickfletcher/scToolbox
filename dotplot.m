function [ax,hs,cb]=dotplot(genes, ncounts, tcounts, glist, ident, ctnames, options)
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
    options.min_prct = 1
    options.prct_leg=3
    options.prct_leg_height=0.05
    options.cb_prct_gap=0.01
    options.cb_gap=0.02
    options.cb_width=0.02
    options.cblabel='mean log_{10} expression'
    options.only_expressing=false
    options.verbose=false
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
[~,EXPR,PRCT]=getExpression(genes(gix,:),ncounts(gix,cellsub),tcounts(gix,cellsub),ident,'only_expressing',options.only_expressing);

[~,ctix]=ismember(ctnames,identnames);
EXPR=EXPR(:,ctix);
PRCT=PRCT(:,ctix);

lowpct=all(PRCT<options.min_prct,2);
EXPR(lowpct,:)=[];
PRCT(lowpct,:)=[];
if options.verbose && nnz(lowpct)>0
    disp("Genes expressed in less than min_prct for all cell groups:")
    disp(G(lowpct)')
end
G(lowpct)=[];
Ggrp(lowpct)=[];

%other filters? ctname-max/min pct list

%other sort filters?

EXPR=log10(EXPR+1);
Glab=strcat('\it',G);
[ax,hs,cb]=plotDotPlot(Glab,ctnames,PRCT,EXPR,fh,Ggrp,options.sortby,options);