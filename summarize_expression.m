function expr_summary = summarize_expression(genes,ncounts,group,options)
arguments
    genes
    ncounts
    group = []
    options.block = []
    options.method = 'mean'
    options.trim_prct = 0
    options.doPrct = true
    options.doPooled = false %true: Cat1, nonCat1, Cat2, ... 
    options.only_expressing = false;
    options.pairwise_contrasts=false
end
%streamlined version of expression summary table (no more thresholds)

%groupsummary????

%remove cats first?
if isempty(group)
    group=ones(size(ncounts,2),1);
end
if ~iscategorical(group)
    group=categorical(group);
end
group=removecats(group);