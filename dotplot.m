function [ax,hs,cb]=dotplot(genes, ncounts, tcounts, glist, ident, ctnames, options)
arguments
    genes
    ncounts
    tcounts
    glist
    ident
    ctnames
    options.Ggroup=[]
    options.sortby='maxsize'
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

ident=removecats(ident);
identnames=categories(ident);

[gix, G]=getGeneIndices(glist,genes.name);
[~,EXPR,PRCT]=getExpression(genes(gix,:),ncounts(gix,:),tcounts(gix,:),ident,'only_expressing',options.only_expressing);

EXPR=EXPR(:,ismember(identnames,ctnames));
PRCT=PRCT(:,ismember(identnames,ctnames));

EXPR=log10(EXPR+1);

Glab=strcat('\it',G);

[ax,hs,cb]=plotDotPlot(Glab,ctnames,PRCT,EXPR,fh,options.Ggroup,options.sortby,options);